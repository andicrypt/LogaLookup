use crate::{PolyIOP, PolyIOPErrors, PolynomialCommitmentScheme, SumCheck, ZeroCheck};
use arithmetic::VPAuxInfo;
use ark_ec::pairing::Pairing;
use ark_ff::{One, PrimeField, Zero};
use ark_poly::DenseMultilinearExtension;
use dashmap::DashMap;
use std::sync::Arc;
use transcript::IOPTranscript;

mod util;
use self::util::*;
use super::structs::LogaPreprocessedTable;

pub trait LookupCheck<E, PCS>: ZeroCheck<E::ScalarField>
where
    E: Pairing,
    PCS: PolynomialCommitmentScheme<E>,
{
    type LookupCheckSubClaim;
    type LookupCheckProof;

    fn init_transcript() -> Self::Transcript;

    fn preprocess_table(
        pcs_param: &PCS::ProverParam,
        nv: usize,
        t: &[E::ScalarField],
    ) -> Result<LogaPreprocessedTable<E, PCS>, PolyIOPErrors>;

    #[allow(clippy::type_complexity)]
    fn prove(
        pcs_param: &PCS::ProverParam,
        f: &Self::MultilinearExtension,
        preprocessed_table: &LogaPreprocessedTable<E, PCS>,
        transcript: &mut IOPTranscript<E::ScalarField>,
    ) -> Result<
        (
            Self::LookupCheckProof,
            Self::MultilinearExtension,
            Self::MultilinearExtension,
            Self::MultilinearExtension,
            Self::MultilinearExtension,
        ),
        PolyIOPErrors,
    >;

    fn verify(
        proof: &Self::LookupCheckProof,
        // zc_aux_info: &VPAuxInfo<E::ScalarField>,
        sc_aux_info: &VPAuxInfo<E::ScalarField>,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::LookupCheckSubClaim, PolyIOPErrors>;
}

/// A lookup check subclaim consists of
/// - A zero check subclaim
/// - A value beta
pub struct LookupCheckSubClaim<F: PrimeField, ZC: ZeroCheck<F>> {
    // pub zero_check_subclaim: ZC::ZeroCheckSubClaim,
    pub sum_subclaim: <ZC as SumCheck<F>>::SumCheckSubClaim,

    /// Challenges beta and alpha
    pub challenges: (F, F),

    // the challenge r which is used to build eq(x, r)
    pub rchallenges: Vec<F>,
}

pub struct LookupCheckProof<E, PCS, ZC>
where
    E: Pairing,
    PCS: PolynomialCommitmentScheme<E>,
    ZC: ZeroCheck<E::ScalarField>,
{
    // pub zc_proof: ZC::ZeroCheckProof,
    pub sc_proof: <ZC as SumCheck<E::ScalarField>>::SumCheckProof,
    pub f_comm: PCS::Commitment,
    pub m_comm: PCS::Commitment,
    pub a_comm: PCS::Commitment,
    pub b_comm: PCS::Commitment,
    pub h_comm: PCS::Commitment,
}

impl<E, PCS> LookupCheck<E, PCS> for PolyIOP<E::ScalarField>
where
    E: Pairing,
    PCS: PolynomialCommitmentScheme<E, Polynomial = Arc<DenseMultilinearExtension<E::ScalarField>>>,
{
    type LookupCheckSubClaim = LookupCheckSubClaim<E::ScalarField, Self>;
    type LookupCheckProof = LookupCheckProof<E, PCS, Self>;

    fn init_transcript() -> Self::Transcript {
        IOPTranscript::<E::ScalarField>::new(b"Initializing LookupCheck transcript")
    }

    fn preprocess_table(
        pcs_param: &PCS::ProverParam,
        nv: usize,
        table: &[E::ScalarField],
    ) -> Result<LogaPreprocessedTable<E, PCS>, PolyIOPErrors> {
        if table.len() != 1 << nv {
            return Err(PolyIOPErrors::InvalidParameters(
                "table does not have sufficient elements".to_string(),
            ));
        }
        let t = Arc::new(
            DenseMultilinearExtension::<E::ScalarField>::from_evaluations_slice(nv, table),
        );
        let t_comm = PCS::commit(pcs_param, &t)?;
        let h_t = DashMap::<E::ScalarField, E::ScalarField>::new();
        for val in table.iter() {
            *h_t.entry(*val).or_insert_with(E::ScalarField::zero) += E::ScalarField::one();
        }
        Ok(LogaPreprocessedTable {
            table: Vec::from(table),
            table_map: h_t,
            t,
            t_comm,
        })
    }

    fn prove(
        pcs_param: &PCS::ProverParam,
        f: &Self::MultilinearExtension,
        preprocessed_table: &LogaPreprocessedTable<E, PCS>,
        transcript: &mut IOPTranscript<E::ScalarField>,
    ) -> Result<
        (
            Self::LookupCheckProof,
            Self::MultilinearExtension,
            Self::MultilinearExtension,
            Self::MultilinearExtension,
            Self::MultilinearExtension,
        ),
        PolyIOPErrors,
    > {
        // Commit to f(x)
        let f_comm = PCS::commit(pcs_param, f)?;
        transcript.append_serializable_element(b"f_comm", &f_comm)?;

        // Compute and commit to m(x)
        let m_poly =
            compute_multiplicity_poly(f, &preprocessed_table.t, &preprocessed_table.table_map)?;
        let m_comm = PCS::commit(pcs_param, &m_poly)?;
        transcript.append_serializable_element(b"m_comm", &m_comm)?;

        // Challenge phase
        let beta = transcript.get_and_append_challenge(b"beta")?;

        // Compute and commit A(x), B(x) such that
        //      A(x) = m(x) / (beta + t(x))
        //      B(x) = 1 / (beta + f(x))
        let a_poly = compute_a(&m_poly, &preprocessed_table.t, &beta)?;
        let b_poly = compute_b(f, &beta)?;

        let a_comm = PCS::commit(pcs_param, &a_poly)?;
        let b_comm = PCS::commit(pcs_param, &b_poly)?;

        transcript.append_serializable_element(b"a_comm", &a_comm)?;
        transcript.append_serializable_element(b"b_comm", &b_comm)?;

        // Compute H as polynomial of mu+1 variables such that
        //      H(X,x) := (1-X).p(x) + X.q(x)
        //      where p(x) = A(x).(beta + t(x)) - m(x)
        //            q(x) = B(x).(beta + f(x)) - 1
        let h_poly = compute_h(&a_poly, &b_poly, &f, &preprocessed_table.t, &m_poly, &beta)?;
        let h_comm = PCS::commit(pcs_param, &h_poly)?;

        transcript.append_serializable_element(b"h_comm", &h_comm)?;

        let length = h_poly.num_vars;
        let r = transcript.get_and_append_challenge_vectors(b"0check r", length)?;

        let h_hat = build_h_hat(&h_poly, r.as_ref())?;


        // Compute L as A(x) - B(x) and
        // lift the domain of L to mu + 1
        let l_poly = build_l_virtual_lifted(&a_poly, &b_poly)?;

        let gamma = transcript.get_and_append_challenge(b"gamma")?;

        let batch_poly = &h_hat + &(l_poly * gamma);

        let sc_proof = <Self as SumCheck<E::ScalarField>>::prove(&batch_poly, transcript)?;

        Ok((
            LookupCheckProof {
                sc_proof,
                f_comm,
                m_comm,
                a_comm,
                b_comm,
                h_comm,
            },
            m_poly,
            a_poly,
            b_poly,
            h_poly,
        ))
    }

    fn verify(
        proof: &Self::LookupCheckProof,
        sc_aux_info: &VPAuxInfo<E::ScalarField>,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::LookupCheckSubClaim, PolyIOPErrors> {
        // update transcript and generate challenge
        transcript.append_serializable_element(b"f_comm", &proof.f_comm)?;
        transcript.append_serializable_element(b"m_comm", &proof.m_comm)?;
        let beta = transcript.get_and_append_challenge(b"beta")?;

        transcript.append_serializable_element(b"a_comm", &proof.a_comm)?;
        transcript.append_serializable_element(b"b_comm", &proof.b_comm)?;
        transcript.append_serializable_element(b"h_comm", &proof.h_comm)?;

        let length = sc_aux_info.num_variables;
        let r = transcript.get_and_append_challenge_vectors(b"0check r", length)?;

        let gamma: <E as Pairing>::ScalarField = transcript.get_and_append_challenge(b"gamma")?;

        let sum_subclaim = <Self as SumCheck<E::ScalarField>>::verify(
            E::ScalarField::zero(),
            &proof.sc_proof,
            &sc_aux_info,
            transcript,
        )?;

        Ok(LookupCheckSubClaim {
            sum_subclaim,
            challenges: (beta, gamma),
            rchallenges: r,
        })
    }
}

#[cfg(test)]
mod test {
    use super::LookupCheck;
    use super::LookupCheckSubClaim;
    use crate::poly_iop::lookup::structs::LogaPreprocessedTable;
    use crate::{
        pcs::{prelude::MultilinearKzgPCS, PolynomialCommitmentScheme},
        poly_iop::{errors::PolyIOPErrors, PolyIOP},
    };
    use arithmetic::{eq_eval, VPAuxInfo};
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_ec::pairing::Pairing;
    use ark_ff::PrimeField;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::rand::Rng;
    use ark_std::test_rng;
    use ark_std::UniformRand;

    use std::marker::PhantomData;
    use std::sync::Arc;

    // Verify A(x)⋅(β+t(x))=m(x).
    // Verify B(x)⋅(β+f(x))=1.
    fn check_a_b_on_boolean_point<F: PrimeField>(
        a_poly: &Arc<DenseMultilinearExtension<F>>,
        b_poly: &Arc<DenseMultilinearExtension<F>>,
        t_poly: &Arc<DenseMultilinearExtension<F>>,
        f_poly: &Arc<DenseMultilinearExtension<F>>,
        m_poly: &Arc<DenseMultilinearExtension<F>>,
        beta: F,
        num_samples: usize,
    ) {
        let mut rng = test_rng();
        let num_vars = a_poly.num_vars;
        assert_eq!(num_vars, b_poly.num_vars);
        assert_eq!(num_vars, t_poly.num_vars);
        assert_eq!(num_vars, f_poly.num_vars);
        assert_eq!(num_vars, m_poly.num_vars);

        for _ in 0..num_samples {
            let point: Vec<F> = (0..num_vars)
                .map(|_| {
                    if bool::rand(&mut rng) {
                        F::one()
                    } else {
                        F::zero()
                    }
                })
                .collect();

            // Evaluate all polynomials at that Boolean point
            let a_eval = a_poly.evaluate(&point).unwrap();
            let b_eval = b_poly.evaluate(&point).unwrap();
            let t_eval = t_poly.evaluate(&point).unwrap();
            let f_eval = f_poly.evaluate(&point).unwrap();
            let m_eval = m_poly.evaluate(&point).unwrap();

            assert_eq!(
                a_eval * (beta + t_eval),
                m_eval,
                "a_eval * (beta + f_eval) == m_eval"
            );
            assert_eq!(
                b_eval * (beta + f_eval),
                F::one(),
                "b_eval * (beta + f_eval) == F::one()"
            );
        }
    }

    // Verify p(x)=A(x)⋅(β+t(x))−m(x)=0.
    // Verify q(x)=B(x)⋅(β+f(x))−1=0.
    fn check_p_q_on_boolean_point<F: PrimeField>(
        a_poly: &Arc<DenseMultilinearExtension<F>>,
        b_poly: &Arc<DenseMultilinearExtension<F>>,
        t_poly: &Arc<DenseMultilinearExtension<F>>,
        f_poly: &Arc<DenseMultilinearExtension<F>>,
        m_poly: &Arc<DenseMultilinearExtension<F>>,
        beta: F,
        num_samples: usize,
    ) -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();
        let num_vars = a_poly.num_vars;
        assert_eq!(num_vars, b_poly.num_vars);
        assert_eq!(num_vars, t_poly.num_vars);
        assert_eq!(num_vars, f_poly.num_vars);
        assert_eq!(num_vars, m_poly.num_vars);

        for _ in 0..num_samples {
            // Sample a random point in {0,1}^num_vars
            let point: Vec<F> = (0..num_vars)
                .map(|_| {
                    if bool::rand(&mut rng) {
                        F::one()
                    } else {
                        F::zero()
                    }
                })
                .collect();

            // Evaluate all polynomials at that Boolean point
            let a_eval = a_poly.evaluate(&point).unwrap();
            let b_eval = b_poly.evaluate(&point).unwrap();
            let t_eval = t_poly.evaluate(&point).unwrap();
            let f_eval = f_poly.evaluate(&point).unwrap();
            let m_eval = m_poly.evaluate(&point).unwrap();

            // Compute p(x) = A(x)(β+t(x)) - m(x)
            let p_eval = a_eval * (beta + t_eval) - m_eval;
            if !p_eval.is_zero() {
                return Err(PolyIOPErrors::InvalidProof(
                    "p(x) != 0 at Boolean point".to_string(),
                ));
            }

            // Compute q(x) = B(x)(β+f(x)) - 1
            let q_eval = b_eval * (beta + f_eval) - F::one();
            if !q_eval.is_zero() {
                return Err(PolyIOPErrors::InvalidProof(
                    "q(x) != 0 at Boolean point".to_string(),
                ));
            }
        }

        Ok(())
    }

    fn test_lookup_check_helper<E, PCS>(
        f: &Arc<DenseMultilinearExtension<E::ScalarField>>,
        preprocessed_table: &LogaPreprocessedTable<E, PCS>,
        pcs_param: &PCS::ProverParam,
    ) -> Result<(), PolyIOPErrors>
    where
        E: Pairing,
        PCS: PolynomialCommitmentScheme<
            E,
            Polynomial = Arc<DenseMultilinearExtension<E::ScalarField>>,
        >,
    {
        // 1) Generate proof
        let mut transcript = <PolyIOP<E::ScalarField> as LookupCheck<E, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;

        let (proof, m_poly, a_poly, b_poly, h_poly) =
            <PolyIOP<E::ScalarField> as LookupCheck<E, PCS>>::prove(
                pcs_param,
                f,
                preprocessed_table,
                &mut transcript,
            )?;

        // 2) Verify the proof
        let mut transcript = <PolyIOP<E::ScalarField> as LookupCheck<E, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;

        // max_degree = 2, num_variables = f.num_vars() + 1 for H, L in F_{mu+1}
        let sc_aux_info: VPAuxInfo<E::ScalarField> = VPAuxInfo {
            max_degree: 2,
            num_variables: f.num_vars() + 1,
            phantom: PhantomData::default(),
        };

        let LookupCheckSubClaim {
            sum_subclaim,
            challenges,
            rchallenges,
        } = <PolyIOP<E::ScalarField> as LookupCheck<E, PCS>>::verify(
            &proof,
            &sc_aux_info,
            &mut transcript,
        )?;

        // 3) Check A and B relations
        check_a_b_on_boolean_point(
            &a_poly,
            &b_poly,
            &preprocessed_table.t,
            f,
            &m_poly,
            challenges.0, // beta
            10,
        );

        // 4) Check H and L relations
        check_p_q_on_boolean_point(
            &a_poly,
            &b_poly,
            &preprocessed_table.t,
            f,
            &m_poly,
            challenges.0, // beta
            10,
        )?;

        let h_eval = h_poly
            .evaluate(&sum_subclaim.point[..h_poly.num_vars()])
            .unwrap();
        let eq_x_r_eval = eq_eval(&sum_subclaim.point, &rchallenges)?;
        let h_hat_value = h_eval * eq_x_r_eval;

        // 5) Check sum-check subclaim
        let a_eval = a_poly
            .evaluate(&sum_subclaim.point[..a_poly.num_vars()])
            .unwrap();
        let b_eval = b_poly
            .evaluate(&sum_subclaim.point[..b_poly.num_vars()])
            .unwrap();
        let l_eval = a_eval - b_eval;
        let batch_poly_eval = h_hat_value + challenges.1 * l_eval;
        assert_eq!(
            batch_poly_eval, sum_subclaim.expected_evaluation,
            "sumcheck on batch_poly not satisfied"
        );

        Ok(())
    }

    fn test_lookup_check(num_vars: usize) -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        // 1) Generate srs
        let srs = MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, num_vars + 1)?;
        let (pcs_param, _) = MultilinearKzgPCS::<Bls12_381>::trim(srs, None, Some(num_vars + 1))?;

        // 2) Generate the table, where the last half is padded.
        let half_n = 1 << (num_vars - 1);
        let half_table = DenseMultilinearExtension::<Fr>::rand(num_vars - 1, &mut rng);

        let mut table = half_table.evaluations;
        table.append(&mut vec![table[half_n - 1]; half_n]);

        let preprocessed_table = <PolyIOP<Fr> as LookupCheck<
            Bls12_381,
            MultilinearKzgPCS<Bls12_381>,
        >>::preprocess_table(&pcs_param, num_vars, &table)?;

        // 3) Generate lookups
        let lookups = (0..half_n)
            .map(|i| vec![table[i]; 2])
            .collect::<Vec<_>>()
            .concat();
        let f = Arc::new(DenseMultilinearExtension::<Fr>::from_evaluations_vec(
            num_vars, lookups,
        ));

        // 4) Generate proof & verify proof
        test_lookup_check_helper::<Bls12_381, MultilinearKzgPCS<Bls12_381>>(
            &f,
            &preprocessed_table,
            &pcs_param,
        )?;

        Ok(())
    }

    fn test_lookup_check_bad_path(num_vars: usize) -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        // 1) Generate SRS
        let srs =
            MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, num_vars + 1).unwrap();
        let (pcs_param, _) =
            MultilinearKzgPCS::<Bls12_381>::trim(srs, None, Some(num_vars + 1)).unwrap();

        // 2) Generate lookup table
        let half_n = 1 << (num_vars - 1);
        let half_table = DenseMultilinearExtension::<Fr>::rand(num_vars - 1, &mut rng);
        let mut table = half_table.evaluations;
        table.append(&mut vec![table[half_n - 1]; half_n]);

        let preprocessed_table = <PolyIOP<Fr> as LookupCheck<
            Bls12_381,
            MultilinearKzgPCS<Bls12_381>,
        >>::preprocess_table(&pcs_param, num_vars, &table)
        .unwrap();

        // 3) Generate correct lookups
        let mut lookups = (0..half_n)
            .map(|i| vec![table[i]; 2])
            .collect::<Vec<_>>()
            .concat();

        // 4) Inject a random mistake into f
        let tamper_index = rng.gen_range(0..lookups.len());
        lookups[tamper_index] += Fr::from(5u64); // small random corruption

        let f = Arc::new(DenseMultilinearExtension::<Fr>::from_evaluations_vec(
            num_vars, lookups,
        ));

        // 5) Generate proof & try to verify — expected to fail
        let result = test_lookup_check_helper::<Bls12_381, MultilinearKzgPCS<Bls12_381>>(
            &f,
            &preprocessed_table,
            &pcs_param,
        );

        match result {
            Ok(_) => panic!("Verification should fail because f(x) was tampered"),
            Err(e) => {
                println!("Verification failed as expected. Error: {:?}", e);
            },
        }

        Ok(())
    }

    #[test]
    fn test_1() -> Result<(), PolyIOPErrors> {
        test_lookup_check(1)?;
        test_lookup_check_bad_path(1)?;
        Ok(())
    }

    #[test]
    fn test_10() -> Result<(), PolyIOPErrors> {
        test_lookup_check(10)?;
        test_lookup_check_bad_path(10)?;
        Ok(())
    }

    #[test]
    fn test_15() -> Result<(), PolyIOPErrors> {
        test_lookup_check(15)?;
        test_lookup_check_bad_path(15)?;
        Ok(())
    }
}
