#![allow(warnings)]
pub mod field;
pub mod ntt;

use std::arch::x86_64::*;
use crate::field::fields::*;
use crate::ntt::ntt::*;
use crate::ntt::intt::*;
use crate::ntt::transpose::*;

fn main() {
    println!("Hello, world!");
}

// #[cfg(test)]
mod tests {
    use std::arch::x86_64::*;
    use crate::field::fields::*;
    use crate::ntt::ntt::*;
    use crate::ntt::intt::*;
    use crate::ntt::transpose::*;

    unsafe fn load_i16_to_m256(input: &[i16; 256]) -> [__m256i; 16] {
        let mut output: [__m256i; 16] = [_mm256_setzero_si256(); 16];
        for i in 0..16 {
            output[i] = _mm256_loadu_si256(input.as_ptr().add(i * 16) as *const __m256i);
        }
        output
    }

    macro_rules! test_ntt_roundtrip {
        ($func_name:ident, $q:expr, $ntt_fn:ident, $intt_fn:ident) => {
            #[test]
            fn $func_name() {
                if !is_x86_feature_detected!("avx2") { return; }
                
                unsafe {
                    let q = $q as i32;
                    // 注意：Scaling factor N 通常是多項式長度 256
                    let scaling_factor = 32i32;

                    // 1. 建立測資
                    let mut p_data = [0i16; 256];
                    for i in 0..256 {
                        p_data[i] = (1 as i32 % q) as i16;
                    }
                    let p_original = p_data.clone();

                    // 2. 載入至向量 (使用 read_unaligned 處理一般 i16 陣列)
                    let ptr = p_data.as_ptr() as *const [__m256i; 16];
                    let input_vecs = std::ptr::read_unaligned(ptr);

                    // 3. 執行正向與逆向變換
                    let ntt_res = crate::ntt::ntt::$ntt_fn(input_vecs);
                    let intt_res = crate::ntt::intt::$intt_fn(ntt_res);

                    // 4. 寫回陣列
                    let mut p_actual = [0i16; 256];
                    let out_ptr = p_actual.as_mut_ptr() as *mut [__m256i; 16];
                    std::ptr::write_unaligned(out_ptr, intt_res);

                    // 5. 驗證同餘關係
                    for i in 0..256 {
                        let a = p_original[i] as i32;
                        let b = p_actual[i] as i32;
                        
                        // 理論值: (a * scaling_factor) mod q
                        let expected_mod_q = (a * scaling_factor).rem_euclid(q);
                        // 實際值取模: b mod q
                        let actual_mod_q = b.rem_euclid(q);

                        // 判斷是否同餘: 只要 actual % q == expected % q 即可
                        // if actual_mod_q != expected_mod_q {
                        //     println!("\n=== DEBUG: Round-trip FAILED (Congruence Check) for q={} ===", q);
                        //     println!("Index {} failed.", i);
                        //     println!("Original[{}]: {}", i, a);
                        //     println!("Actual (Raw) : {}", b);
                        //     println!("Actual (mod q): {}", actual_mod_q);
                        //     println!("Expected (mod q): {}", expected_mod_q);
                            
                        //     println!("\n--- Full Polynomial Arrays ---");
                        //     println!("Original: {:?}\n", p_original);
                        //     println!("Actual Raw: {:?}\n", p_actual);
                            
                        //     let p_expected_mod: Vec<i16> = p_original.iter()
                        //         .map(|&x| (x as i32 * scaling_factor).rem_euclid(q) as i16)
                        //         .collect();
                        //     println!("Expected (mod q): {:?}", p_expected_mod);
                        //     println!("==========================================\n");

                        //     panic!(
                        //         "Round-trip FAILED at index {}: {} is not congruent to {} (mod {})", 
                        //         i, b, expected_mod_q, q
                        //     );
                        // }
                        let mut error_indices = Vec::new();
                        for i in 0..256 {
                            let a = p_original[i] as i32;
                            let b = p_actual[i] as i32;
                            
                            let expected_mod_q = (a * scaling_factor).rem_euclid(q);
                            let actual_mod_q = b.rem_euclid(q);

                            if actual_mod_q != expected_mod_q {
                                error_indices.push(i);
                            }
                        }
                        if !error_indices.is_empty() {
                            println!("\n=== DEBUG: Round-trip FAILED (Congruence Check) for q={} ===", q);
                            println!("Total errors found: {}/256", error_indices.len());
                            println!("{:<10} | {:<10} | {:<10} | {:<10}", "Index", "Original", "Expected", "Actual(Raw)");
                            println!("{:-<50}", "");

                            for &i in &error_indices {
                                let a = p_original[i] as i32;
                                let mut b = p_actual[i] as i32;
                                let expected = (a * scaling_factor).rem_euclid(q);
                                
                                // 僅列印出錯的位元
                                println!("{:<10} | {:<10} | {:<10} | {:<10}", i, a, expected, b);
                            }
                            println!("==========================================\n");

                            panic!(
                                "Round-trip FAILED for q={}: {} mismatches found.", 
                                q, error_indices.len()
                            );
                        }
                    }
                    println!("SUCCESS: Round-trip for q={} passed (Congruent within modulo).", q);
                }
            }
        };
    }
    test_ntt_roundtrip!(roundtrip_7681, 7681, ntt_7681, intt_7681);
    // test_ntt_roundtrip!(roundtrip_10753, 10753, ntt_10753, intt_10753);
    // test_ntt_roundtrip!(roundtrip_11777, 11777, ntt_11777, intt_11777);
    // test_ntt_roundtrip!(roundtrip_12289, 12289, ntt_12289, intt_12289);
    // test_ntt_roundtrip!(roundtrip_13313, 13313, ntt_13313, intt_13313);
    // test_ntt_roundtrip!(roundtrip_15361, 15361, ntt_15361, intt_15361);
}
