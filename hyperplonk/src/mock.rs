// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

use arithmetic::identity_permutation;
use ark_ff::PrimeField;
use ark_std::{log2, test_rng};

use crate::{
    custom_gate::CustomizedGates,
    selectors::SelectorColumn,
    structs::{HyperPlonkIndex, HyperPlonkParams},
    witness::WitnessColumn,
};

#[derive(Debug)]
pub struct MockCircuit<F: PrimeField> {
    pub public_inputs: Vec<F>,
    pub witnesses: Vec<WitnessColumn<F>>,
    pub index: HyperPlonkIndex<F>,
}

impl<F: PrimeField> MockCircuit<F> {
    /// Number of variables in a multilinear system
    pub fn num_variables(&self) -> usize {
        self.index.num_variables()
    }

    /// Number of constraint rows
    pub fn num_constraints(&self) -> usize {
        self.index.params.num_constraints
    }

    /// number of selector columns
    pub fn num_selector_columns(&self) -> usize {
        self.index.num_selector_columns()
    }

    /// number of lookup selector columns
    pub fn num_lk_selector_columns(&self) -> usize {
        self.index.num_lk_selector_columns()
    }

    /// number of witness columns
    pub fn num_witness_columns(&self) -> usize {
        self.index.num_witness_columns()
    }
}

impl<F: PrimeField> MockCircuit<F> {
    /// Generate a mock plonk circuit for the input constraint size.
    pub fn new(
        num_constraints: usize,
        gate: &CustomizedGates,
        lk_gate: &CustomizedGates,
        full_table: bool,
    ) -> MockCircuit<F> {
        let mut rng = test_rng();
        let nv = log2(num_constraints);
        let num_selectors = gate.num_selector_columns();
        let num_lk_selectors = lk_gate.num_selector_columns();
        let num_witnesses = gate.num_witness_columns();
        let log_n_wires = log2(num_witnesses);
        let merged_nv = nv + log_n_wires;

        let mut selectors: Vec<SelectorColumn<F>> = vec![SelectorColumn::default(); num_selectors];
        let mut lk_selectors: Vec<SelectorColumn<F>> =
            vec![SelectorColumn::default(); num_lk_selectors];
        let mut witnesses: Vec<WitnessColumn<F>> = vec![WitnessColumn::default(); num_witnesses];

        let mut table;
        if full_table {
            table = vec![F::default(); num_constraints]
        } else {
            table = vec![F::default(); num_constraints - 1]
        };

        for current_row in 0..num_constraints {
            // 1) Generate witness for this row.
            let cur_witness: Vec<F> = (0..num_witnesses).map(|_| F::rand(&mut rng)).collect();
            for i in 0..num_witnesses {
                witnesses[i].append(cur_witness[i]);
            }

            // 2) Generate selectors for this row such that row evaluates to 0.
            let mut cur_selectors: Vec<F> = (0..(num_selectors - 1))
                .map(|_| F::rand(&mut rng))
                .collect();

            let mut last_selector = F::zero();
            for (index, (coeff, q, wit)) in gate.gates.iter().enumerate() {
                if index != num_selectors - 1 {
                    let mut cur_monomial = if *coeff < 0 {
                        -F::from((-coeff) as u64)
                    } else {
                        F::from(*coeff as u64)
                    };
                    cur_monomial = match q {
                        Some(p) => cur_monomial * cur_selectors[*p],
                        None => cur_monomial,
                    };
                    for wit_index in wit.iter() {
                        cur_monomial *= cur_witness[*wit_index];
                    }
                    last_selector += cur_monomial;
                } else {
                    let mut cur_monomial = if *coeff < 0 {
                        -F::from((-coeff) as u64)
                    } else {
                        F::from(*coeff as u64)
                    };
                    for wit_index in wit.iter() {
                        cur_monomial *= cur_witness[*wit_index];
                    }
                    last_selector /= -cur_monomial;
                }
            }
            cur_selectors.push(last_selector);
            for i in 0..num_selectors {
                selectors[i].append(cur_selectors[i]);
            }

            // 3) Generate lookup selectors, evaluate the row and push to table.

            let mut lookup_value = F::zero();
            let mut cur_lk_selectors: Vec<F>;
            if current_row != num_constraints - 1 {
                cur_lk_selectors = (0..(num_lk_selectors)).map(|_| F::rand(&mut rng)).collect();

                for (coeff, q_lk, wit) in lk_gate.gates.iter() {
                    let mut cur_lk_monomial = if *coeff < 0 {
                        -F::from((-coeff) as u64)
                    } else {
                        F::from(*coeff as u64)
                    };
                    cur_lk_monomial = match q_lk {
                        Some(p) => cur_lk_monomial * cur_lk_selectors[*p],
                        None => cur_lk_monomial,
                    };
                    for wit_index in wit.iter() {
                        cur_lk_monomial *= cur_witness[*wit_index];
                    }
                    lookup_value += cur_lk_monomial;
                }
                table[current_row] = lookup_value;
            } else {
                cur_lk_selectors = (0..(num_lk_selectors - 1))
                    .map(|_| F::rand(&mut rng))
                    .collect();
                let mut last_lk_selector = F::zero();
                for (index, (coeff, q_lk, wit)) in lk_gate.gates.iter().enumerate() {
                    if index != num_lk_selectors - 1 {
                        let mut cur_lk_monomial = if *coeff < 0 {
                            -F::from((-coeff) as u64)
                        } else {
                            F::from(*coeff as u64)
                        };
                        cur_lk_monomial = match q_lk {
                            Some(p) => cur_lk_monomial * cur_lk_selectors[*p],
                            None => cur_lk_monomial,
                        };
                        for wit_index in wit.iter() {
                            cur_lk_monomial *= cur_witness[*wit_index];
                        }
                        last_lk_selector += cur_lk_monomial;
                    } else {
                        let mut cur_lk_monomial = if *coeff < 0 {
                            -F::from((-coeff) as u64)
                        } else {
                            F::from(*coeff as u64)
                        };
                        for wit_index in wit.iter() {
                            cur_lk_monomial *= cur_witness[*wit_index];
                        }

                        // acc_last_selector + x * cur_lk_monomial = table[current - 1]
                        last_lk_selector =
                            (table[current_row - 1] - last_lk_selector) / cur_lk_monomial;
                    }
                }
                cur_lk_selectors.push(last_lk_selector);
            }

            for i in 0..num_lk_selectors {
                lk_selectors[i].append(cur_lk_selectors[i]);
            }
        }
        let pub_input_len = ark_std::cmp::min(4, num_constraints);
        let public_inputs: Vec<F> = witnesses[0].0[0..pub_input_len].to_vec();

        let params = HyperPlonkParams {
            num_constraints,
            num_pub_input: public_inputs.len(),
            gate_func: gate.clone(),
            lk_gate_func: lk_gate.clone(),
        };

        let permutation = identity_permutation(merged_nv as usize, 1);

        let index = HyperPlonkIndex {
            params,
            permutation,
            selectors: selectors.clone(),
            lk_selectors: lk_selectors.clone(),
            table,
        };

        Self {
            public_inputs,
            witnesses,
            index,
        }
    }

    pub fn is_satisfied(&self) -> bool {
        for current_row in 0..self.num_constraints() {
            // gate func satisfaction
            let mut cur = F::zero();

            for (coeff, q, wit) in self.index.params.gate_func.gates.iter() {
                let mut cur_monomial = if *coeff < 0 {
                    -F::from((-coeff) as u64)
                } else {
                    F::from(*coeff as u64)
                };
                cur_monomial = match q {
                    Some(p) => cur_monomial * self.index.selectors[*p].0[current_row],
                    None => cur_monomial,
                };
                for wit_index in wit.iter() {
                    cur_monomial *= self.witnesses[*wit_index].0[current_row];
                }
                cur += cur_monomial;
            }
            if !cur.is_zero() {
                return false;
            }

            // lookup gate func satisfaction
            cur = F::zero();

            for (coeff, q_lk, wit) in self.index.params.lk_gate_func.gates.iter() {
                let mut cur_monomial = if *coeff < 0 {
                    -F::from((-coeff) as u64)
                } else {
                    F::from(*coeff as u64)
                };
                if let Some(p_lk) = q_lk {
                    cur_monomial *= self.index.lk_selectors[*p_lk].0[current_row];
                }
                for w in wit.iter() {
                    cur_monomial *= self.witnesses[*w].0[current_row];
                }
                cur += cur_monomial;
            }

            if !self.index.table.contains(&cur) {
                return false;
            }
        }

        true
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{errors::HyperPlonkErrors, LogaHyperPlonkSNARK};
    use ark_bls12_381::{Bls12_381, Fr};
    use subroutines::{
        pcs::{
            prelude::{MultilinearKzgPCS, MultilinearUniversalParams},
            PolynomialCommitmentScheme,
        },
        poly_iop::PolyIOP,
    };

    const SUPPORTED_SIZE: usize = 7;
    const MIN_NUM_VARS: usize = 3;
    const MAX_NUM_VARS: usize = 7;
    const CUSTOM_DEGREE: [usize; 3] = [1, 2, 4];

    #[test]
    fn test_sat_function() {
        let nv = 2;

        // q0 * x0^2 - x1 = 0
        let gate_func = CustomizedGates {
            gates: vec![(1, Some(0), vec![0; 2]), (-1, None, vec![1])],
        };
        // qlk0 * x0 + qlk1 * x1 in table
        let lk_gate_func = CustomizedGates {
            gates: vec![(1, Some(0), vec![0]), (1, Some(1), vec![1])],
        };

        let params = HyperPlonkParams {
            num_constraints: 1 << nv,
            num_pub_input: 1,
            gate_func,
            lk_gate_func,
        };
        let index = HyperPlonkIndex {
            params,
            permutation: vec![Fr::from(1); 4],
            selectors: vec![SelectorColumn(vec![Fr::from(1); 4])],
            lk_selectors: vec![SelectorColumn(vec![Fr::from(1); 4]); 2],
            table: vec![Fr::from(2), Fr::from(6), Fr::from(12), Fr::from(20)],
        };

        let w1 = WitnessColumn(vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)]);
        let w2 = WitnessColumn(vec![Fr::from(1), Fr::from(4), Fr::from(9), Fr::from(16)]);
        let circuit = MockCircuit {
            public_inputs: vec![Fr::from(1)],
            witnesses: vec![w1, w2],
            index,
        };

        assert!(circuit.is_satisfied());
    }

    #[test]
    fn test_mock_circuit_sat() {
        for i in 1..5 {
            let vanilla_gate = CustomizedGates::vanilla_plonk_gate();
            let circuit = MockCircuit::<Fr>::new(1 << i, &vanilla_gate, &vanilla_gate, true);
            assert!(circuit.is_satisfied());

            let jf_gate = CustomizedGates::jellyfish_turbo_plonk_gate();
            let circuit = MockCircuit::<Fr>::new(1 << i, &jf_gate, &jf_gate, true);
            assert!(circuit.is_satisfied());

            for num_witness in 2..10 {
                for degree in CUSTOM_DEGREE {
                    let mock_gate = CustomizedGates::mock_gate(num_witness, degree);
                    let circuit = MockCircuit::<Fr>::new(1 << i, &mock_gate, &mock_gate, true);
                    assert!(circuit.is_satisfied());
                    let circuit = MockCircuit::<Fr>::new(1 << i, &mock_gate, &mock_gate, false);
                    assert!(circuit.is_satisfied());
                }
            }
        }
    }

    fn test_mock_circuit_zkp_helper(
        nv: usize,
        gate: &CustomizedGates,
        pcs_srs: &MultilinearUniversalParams<Bls12_381>,
    ) -> Result<(), HyperPlonkErrors> {
        let circuit = MockCircuit::<Fr>::new(1 << nv, gate, gate, true);
        assert!(circuit.is_satisfied());

        let index = circuit.index;
        // generate pk and vks
        let (pk, vk) = <PolyIOP<Fr> as LogaHyperPlonkSNARK<
            Bls12_381,
            MultilinearKzgPCS<Bls12_381>,
        >>::preprocess(&index, pcs_srs)?;
        // generate a proof and verify
        let proof = <PolyIOP<Fr> as LogaHyperPlonkSNARK<
            Bls12_381,
            MultilinearKzgPCS<Bls12_381>,
        >>::prove(&pk, &circuit.public_inputs, &circuit.witnesses)?;

        let verify = <PolyIOP<Fr> as LogaHyperPlonkSNARK<
            Bls12_381,
            MultilinearKzgPCS<Bls12_381>,
        >>::verify(&vk, &circuit.public_inputs, &proof)?;
        assert!(verify);
        Ok(())
    }

    #[test]
    fn test_mock_circuit_zkp() -> Result<(), HyperPlonkErrors> {
        let mut rng = test_rng();
        let pcs_srs =
            MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, SUPPORTED_SIZE)?;
        for nv in MIN_NUM_VARS..MAX_NUM_VARS {
            let vanilla_gate = CustomizedGates::vanilla_plonk_gate();
            test_mock_circuit_zkp_helper(nv, &vanilla_gate, &pcs_srs)?;
        }
        for nv in MIN_NUM_VARS..MAX_NUM_VARS {
            let tubro_gate = CustomizedGates::jellyfish_turbo_plonk_gate();
            test_mock_circuit_zkp_helper(nv, &tubro_gate, &pcs_srs)?;
        }
        let nv = 5;
        for num_witness in 2..6 {
            for degree in CUSTOM_DEGREE {
                let mock_gate = CustomizedGates::mock_gate(num_witness, degree);
                test_mock_circuit_zkp_helper(nv, &mock_gate, &pcs_srs)?;
            }
        }

        Ok(())
    }

    #[test]
    fn test_mock_circuit_e2e() -> Result<(), HyperPlonkErrors> {
        let mut rng = test_rng();
        let pcs_srs =
            MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, SUPPORTED_SIZE)?;
        let nv = MAX_NUM_VARS;

        let turboplonk_gate = CustomizedGates::jellyfish_turbo_plonk_gate();
        test_mock_circuit_zkp_helper(nv, &turboplonk_gate, &pcs_srs)?;

        Ok(())
    }

    #[test]
    fn test_mock_long_selector_e2e() -> Result<(), HyperPlonkErrors> {
        let mut rng = test_rng();
        let pcs_srs =
            MultilinearKzgPCS::<Bls12_381>::gen_srs_for_testing(&mut rng, SUPPORTED_SIZE)?;
        let nv = MAX_NUM_VARS;

        let long_selector_gate = CustomizedGates::super_long_selector_gate();
        test_mock_circuit_zkp_helper(nv, &long_selector_gate, &pcs_srs)?;

        Ok(())
    }
}
