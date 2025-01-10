use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use dashmap::{DashMap, DashSet};
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};
use std::sync::Arc;

use crate::poly_iop::errors::PolyIOPErrors;

pub(super) fn compute_h<F: PrimeField>(
    f: &Arc<DenseMultilinearExtension<F>>,
    t: &[F],
    h_t: &DashMap<F, usize>,
) -> Result<Vec<F>, PolyIOPErrors> {
    let nv = f.num_vars;
    assert!(
        (1 << nv) - 1 == t.len(),
        "lookup size and table size are incorrect"
    );

    let h_f = DashMap::new();

    f.evaluations
        .par_iter()
        .map(|num| -> Result<(), PolyIOPErrors> {
            if h_t.get(num).is_none() {
                return Err(PolyIOPErrors::InvalidProof(format!(
                    "Lookup value {num} is not in table"
                )));
            }
            *h_f.entry(*num).or_insert_with(|| 0) += 1;
            Ok(())
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut evaluations: Vec<F> = vec![];
    let table_set = DashSet::new();
    for num in t.iter() {
        if !table_set.contains(num) {
            let h_t_val: usize = *h_t.get(num).unwrap();
            if let Some(h_f_val) = h_f.get(num) {
                evaluations.append(&mut vec![*num; h_t_val + *h_f_val]);
            } else {
                evaluations.append(&mut vec![*num; h_t_val]);
            }
            table_set.insert(*num);
        }
    }

    Ok(evaluations)
}

pub(super) fn get_primitive_polynomial(nv: usize) -> Result<usize, PolyIOPErrors> {
    if !(2..=33).contains(&nv) {
        return Err(PolyIOPErrors::InvalidProof(
            "Only support primitive polynomial whose degree >= 2 and <= 32".to_string(),
        ));
    }

    let primitive_polynomial: Vec<usize> = vec![
        0b111,                               //2
        0b1011,                              //3
        0b10011,                             //4
        0b100101,                            //5
        0b1000011,                           //6
        0b10000011,                          //7
        0b100011101,                         //8
        0b1000010001,                        //9
        0b10000001001,                       //10
        0b100000000101,                      //11
        0b1000001010011,                     //12
        0b10000000011011,                    //13
        0b100000101000011,                   //14
        0b1000000000000011,                  //15
        0b10001000000001011,                 //16
        0b100000000000001001,                //17
        0b1000000000010000001,               //18
        0b10000000000001100011,              //19
        0b100000000000000001001,             //20
        0b1000000000000000000101,            //21
        0b10000000000000000000011,           //22
        0b100000000000000000100001,          //23
        0b1000000000000000000011011,         //24
        0b10000000000000000000001001,        //25
        0b100000000000000000110000011,       //26
        0b1000000000000000000110000011,      //27
        0b10000000000000000000000001001,     //28
        0b100000000000000000000000000101,    //29
        0b1000000100000000000000000000111,   //30
        0b10000000000000000000000000001001,  //31
        0b100000000100000000000000000000111, //32
    ];
    Ok(primitive_polynomial[nv - 2].clone())
}

/// g_mu(b1,...,b_mu) = (b_mu, b1', ... , b_{mu-1}')
///
pub(super) fn next_element(cur_num: usize, nv: usize) -> usize {
    // Reverse current number and exclude the LSB
    // For example, cur_num = 14 (0b1110)
    // Reverse form: 0b0111
    // After exclusion: 0b111
    let reverse_bit = reverse_bits_func(cur_num, nv);
    let reverse_bit_short = reverse_bit & ((1 << (nv - 1)) - 1); // exclude the first bit

    // Exclude the first and last bit, corresponding to the highest and free term.
    // For example, with primitive polynomial which generate GF(2^3): X^3 + X + 1
    // Binary representation: 1011
    // After exclusion: 01
    let s = get_primitive_polynomial(nv).unwrap();
    let s_trim = (s & ((1 << (nv - 1)) - 1)) >> 1;

    let first_bit = (reverse_bit >> (nv - 1)) & 1;

    let next_num = reverse_bits_func(
        ((reverse_bit_short ^ if first_bit == 1 { s_trim } else { 0 }) << 1) + first_bit,
        nv,
    );

    next_num
}

pub fn embed<F: PrimeField>(
    poly_evals: &[F], //2^num_vars - 1
    nv: usize,
) -> Result<Arc<DenseMultilinearExtension<F>>, PolyIOPErrors> {
    assert!(
        poly_evals.len() == (1 << nv) - 1,
        "Embedded evaluations must be in form of pow(2,nv) - 1"
    );

    let mut embedded_poly_evals: Vec<F> = vec![F::zero(); 1 << nv];
    let mut cur_element: usize = 1 << (nv - 1);

    for item in poly_evals.iter() {
        embedded_poly_evals[cur_element] = *item;
        cur_element = next_element(cur_element, nv);
    }

    Ok(Arc::new(DenseMultilinearExtension::from_evaluations_vec(
        nv,
        embedded_poly_evals,
    )))
}

/// poly_delta(X1,X2...,Xnv) = Xnv * poly(1,X1',...,X_{nv-1}')
///             + (1 - Xnv) * poly(0,X1,...,X_{nv-1})
/// where Xi' := 1 - Xi (if i in S), and Xi' := Xi otherwise.
pub(super) fn compute_poly_delta<F: PrimeField>(
    poly: &Arc<DenseMultilinearExtension<F>>,
    nv: usize,
) -> Result<Arc<DenseMultilinearExtension<F>>, PolyIOPErrors> {
    let s = get_primitive_polynomial(nv).unwrap();

    // trim first and last bit of binary representation of primitive polynomial
    let s_trim = (s & ((1 << (nv - 1)) - 1)) >> 1;

    let poly_evals = &poly.evaluations;

    let mut evaluations: Vec<F> = vec![];
    evaluations.push(poly_evals[0]);

    for i in 1..(1 << nv) {
        let reverse_bit = reverse_bits_func(i, nv);
        let reverse_bit_short = reverse_bit & ((1 << (nv - 1)) - 1); // exclude the first bit

        let x_mu = (reverse_bit >> (nv - 1)) & 1;
        let next: usize;

        next = reverse_bits_func(
            ((reverse_bit_short ^ if x_mu == 1 { s_trim } else { 0 }) << 1) + x_mu,
            nv,
        );

        evaluations.push(poly_evals[next]);
    }

    Ok(Arc::new(DenseMultilinearExtension::from_evaluations_vec(
        nv,
        evaluations,
    )))
}

/// reverse_bits_func reverse the binary representation of number, and truncate its size
/// For example, num=6 (0b110), bit_size=3
/// Reverse form: 0b011
pub fn reverse_bits_func(num: usize, bit_size: usize) -> usize {
    let reversed = num.reverse_bits(); // Reverse all 64 bits
    reversed >> (64 - bit_size) // Shift to keep only the first `bit_size` bits
}

#[cfg(test)]
mod test {
    use super::embed;
    use super::*;
    use ark_bls12_381::Fr;
    use ark_poly::MultilinearExtension;
    use ark_std::test_rng;

    #[test]
    fn test_embed() -> Result<(), PolyIOPErrors> {
        let nv = 3;
        let poly_evals = (1..8).map(Fr::from).collect::<Vec<_>>();

        let poly_embed = embed(&poly_evals, nv)?;

        let poly_embed_evals = &poly_embed.evaluations;

        for x in poly_evals.iter() {
            if !poly_embed_evals.contains(x) {
                return Err(PolyIOPErrors::InvalidParameters(
                    "wrong embedding".to_string(),
                ));
            }
        }
        Ok(())
    }

    #[test]
    fn test_compute_poly_delta() -> Result<(), PolyIOPErrors> {
        let mut rng = test_rng();

        let nv = 3;
        let poly = Arc::new(DenseMultilinearExtension::<Fr>::rand(nv, &mut rng));
        let poly_evals = &poly.evaluations;

        let poly_delta = compute_poly_delta(&poly, nv)?;
        let poly_delta_evals = &poly_delta.evaluations;

        if poly_evals.len() != poly_delta_evals.len() {
            return Err(PolyIOPErrors::InvalidParameters(
                "wrong embedding: Result poly not have the same size".to_string(),
            ));
        }
        for x in poly_evals {
            if !poly_delta_evals.contains(x) {
                return Err(PolyIOPErrors::InvalidParameters(
                    "wrong embedding: Result poly not equal to Original poly".to_string(),
                ));
            }
        }
        Ok(())
    }

    #[test]
    fn test_compute_h() -> Result<(), PolyIOPErrors> {
        let nv = 3;

        // generate the table, whose each element is distinct
        let table = (1..(1 << nv)).map(Fr::from).collect::<Vec<_>>();
        let h_t = DashMap::new();
        for num in table.iter() {
            *h_t.entry(*num).or_insert_with(|| 0) += 1;
        }

        let half_nv = 1 << (nv - 1);
        let lookups = (0..half_nv)
            .map(|i| vec![table[i]; 2])
            .collect::<Vec<_>>()
            .concat();

        let f = Arc::new(DenseMultilinearExtension::from_evaluations_vec(
            nv,
            lookups.clone(),
        ));

        let h = compute_h(&f, &table, &h_t)?;
        if h.len() != lookups.len() + table.len() {
            return Err(PolyIOPErrors::InvalidParameters(
                "wrong computing h: Incorrect result vector size".to_string(),
            ));
        }

        let mut lookups_len_expected = 0;
        let mut cur_table_value = table[0] - Fr::from(1);
        for x in h {
            if x == cur_table_value {
                lookups_len_expected += 1;
            } else {
                cur_table_value = x;
            }
        }

        if lookups_len_expected != lookups.len() {
            return Err(PolyIOPErrors::InvalidParameters(
                "wrong computing h: Incorrect result vector size".to_string(),
            ));
        }

        Ok(())
    }
}
