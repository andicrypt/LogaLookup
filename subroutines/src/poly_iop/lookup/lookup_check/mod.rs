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
            Self::VirtualPolynomial,
            Self::VirtualPolynomial,
        ),
        PolyIOPErrors,
    >;

    fn verify(
        proof: &Self::LookupCheckProof,
        sc_aux_info: &VPAuxInfo<E::ScalarField>,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::LookupCheckSubClaim, PolyIOPErrors>;
}

/// A lookup check subclaim consists of
/// - A zero check subclaim
/// - A value beta
pub struct LookupCheckSubClaim<F: PrimeField, ZC: ZeroCheck<F>> {
    /// Sumcheck subclaim
    pub sum_check_subclaim: <ZC as SumCheck<F>>::SumCheckSubClaim,

    /// Challenges beta, alpha, gamma
    pub challenges: (F, F, F),

    /// r challenge for zerocheck identity checking
    pub rchallenges: Vec<F>,
}

pub struct LookupCheckProof<E, PCS, ZC>
where
    E: Pairing,
    PCS: PolynomialCommitmentScheme<E>,
    ZC: ZeroCheck<E::ScalarField>,
{
    pub sc_proof: <ZC as SumCheck<E::ScalarField>>::SumCheckProof,
    pub f_comm: PCS::Commitment,
    pub m_comm: PCS::Commitment,
    pub a_comm: PCS::Commitment,
    pub b_comm: PCS::Commitment,
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
            Self::VirtualPolynomial,
            Self::VirtualPolynomial,
        ),
        PolyIOPErrors,
    > {
        // Commit lookup & multiplicity polynomials.
        let f_comm = PCS::commit(pcs_param, f)?;
        transcript.append_serializable_element(b"f_comm", &f_comm)?;

        let m_poly =
            compute_multiplicity_poly(f, &preprocessed_table.t, &preprocessed_table.table_map)?;
        let m_comm = PCS::commit(pcs_param, &m_poly)?;
        transcript.append_serializable_element(b"m_comm", &m_comm)?;

        // Create and commit polynomials A, B
        //      A(x) = m(x) / (beta + t(x))
        //      B(x) = 1 / (beta + f(x))
        let beta = transcript.get_and_append_challenge(b"beta")?;

        let a_poly = compute_a(&m_poly, &preprocessed_table.t, &beta)?;
        let b_poly = compute_b(f, &beta)?;

        let a_comm = PCS::commit(pcs_param, &a_poly)?;
        let b_comm = PCS::commit(pcs_param, &b_poly)?;

        transcript.append_serializable_element(b"a_comm", &a_comm)?;
        transcript.append_serializable_element(b"b_comm", &b_comm)?;

        // Build batched virtual polynomial p + alpha * q
        let alpha = transcript.get_and_append_challenge(b"alpha")?;
        let h_poly_virtual = build_h_virtual(
            &a_poly,
            &b_poly,
            f,
            &preprocessed_table.t,
            &m_poly,
            &alpha,
            &beta,
        )?;

        // Build H_hat(X) = M(X) * eq_r(X)
        // where eq_r(X) = r_i * X_i + (1 - r_i) * (1 - X_i)
        let r = transcript.get_and_append_challenge_vectors(b"0check r", f.num_vars)?;
        let h_hat_poly_virtual = h_poly_virtual.build_f_hat(r.as_ref())?;

        // Build L(X) = A(X) - B(X) over boolean hypercube
        let l_poly_virtual = build_l_virtual(&a_poly, &b_poly)?;

        // Batch two sumcheck on H_hat and L by using
        // randomly sampled challenge gamma from verifier
        let gamma = transcript.get_and_append_challenge(b"gamma")?;
        let t_poly_virtual = &h_hat_poly_virtual + &(&l_poly_virtual * gamma);

        let sc_proof = <Self as SumCheck<E::ScalarField>>::prove(&t_poly_virtual, transcript)?;

        Ok((
            LookupCheckProof {
                sc_proof,
                f_comm,
                m_comm,
                a_comm,
                b_comm,
            },
            m_poly,
            a_poly,
            b_poly,
            h_poly_virtual,
            l_poly_virtual,
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

        let alpha = transcript.get_and_append_challenge(b"alpha")?;

        let length = sc_aux_info.num_variables;
        let r = transcript.get_and_append_challenge_vectors(b"0check r", length)?;

        let gamma = transcript.get_and_append_challenge(b"gamma")?;

        let sc_sub_claim = <Self as SumCheck<E::ScalarField>>::verify(
            E::ScalarField::zero(),
            &proof.sc_proof,
            sc_aux_info,
            transcript,
        )?;

        Ok(LookupCheckSubClaim {
            sum_check_subclaim: sc_sub_claim,
            challenges: (beta, alpha, gamma),
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
    use ark_ff::One;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::rand::Rng;
    use ark_std::test_rng;

    use std::marker::PhantomData;
    use std::sync::Arc;

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

        let (proof, m_poly, a_poly, b_poly, h_poly_virtual, l_poly_virtual) =
            <PolyIOP<E::ScalarField> as LookupCheck<E, PCS>>::prove(
                pcs_param,
                f,
                preprocessed_table,
                &mut transcript,
            )?;

        // 2) Verify the proof
        let mut transcript = <PolyIOP<E::ScalarField> as LookupCheck<E, PCS>>::init_transcript();
        transcript.append_message(b"testing", b"initializing transcript for testing")?;

        // max_degree = 3, num_var = mu
        let sc_aux_info: VPAuxInfo<E::ScalarField> = VPAuxInfo {
            max_degree: 3,
            num_variables: f.num_vars(),
            phantom: PhantomData::default(),
        };

        let LookupCheckSubClaim {
            sum_check_subclaim: sum_subclaim,
            challenges,
            rchallenges,
        } = <PolyIOP<E::ScalarField> as LookupCheck<E, PCS>>::verify(
            &proof,
            &sc_aux_info,
            &mut transcript,
        )?;

        let a_eval = a_poly.evaluate(&sum_subclaim.point).unwrap();
        let b_eval = b_poly.evaluate(&sum_subclaim.point).unwrap();
        let m_eval = m_poly.evaluate(&sum_subclaim.point).unwrap();
        let t_eval = preprocessed_table.t.evaluate(&sum_subclaim.point).unwrap();
        let f_eval = f.evaluate(&sum_subclaim.point).unwrap();
        let h_eval = h_poly_virtual.evaluate(&sum_subclaim.point).unwrap();
        let l_eval = l_poly_virtual.evaluate(&sum_subclaim.point).unwrap();
        let (beta, alpha, gamma) = challenges;

        let eq_x_r_eval = eq_eval(&sum_subclaim.point, &rchallenges)?;
        let h_hat_eval = h_eval * eq_x_r_eval;

        assert_eq!(
            h_hat_eval + gamma * l_eval,
            sum_subclaim.expected_evaluation,
            "h_hat_eval + gamma * (a_eval - b_eval) == sum_subclaim.expected_evaluation"
        );

        assert_eq!(
            h_eval,
            a_eval * (beta + t_eval) - m_eval + alpha * (b_eval * (beta + f_eval) - E::ScalarField::one()),
            "h_eval = a_eval * (beta + t_eval) - m_eval + alpha * (b_eval * (beta + f_eval) - F::one())"
        );
        Ok(())
    }

    fn test_lookup_check(num_vars: usize) -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        // 1) Generate srs
        let srs = MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, num_vars)?;
        let (pcs_param, _) = MultilinearKzgPCS::<Bls12_381>::trim(srs, None, Some(num_vars))?;

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
