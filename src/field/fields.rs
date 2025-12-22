#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::field::macros::*;
// RNS: [257, 3329, 7681, 7937, 9473, 10753]

// reduce 2^32-99
#[target_feature(enable = "avx2")]
pub unsafe fn reduce32(even_64: __m256i, odd_64: __m256i) -> __m256i {
    
    let reduce_lane = |x: __m256i| -> __m256i {
        let c99 = _mm256_set1_epi64x(99);
        let p_val = _mm256_set1_epi64x(4294967197); // 2^32 - 99
        let mask_low = _mm256_set1_epi64x(0x00000000FFFFFFFF);
        let zero = _mm256_setzero_si256();

        // Round 1
        let h1 = _mm256_srai_epi64::<32>(x); 
        let l1 = _mm256_and_si256(x, mask_low);
        let sum1 = _mm256_add_epi64(_mm256_mul_epi32(h1, c99), l1);

        // Round 2 (Convergence)
        let h2 = _mm256_srai_epi64::<32>(sum1);
        let l2 = _mm256_and_si256(sum1, mask_low);
        let sum2 = _mm256_add_epi64(_mm256_mul_epi32(h2, c99), l2);

        // Correction: if sum2 < 0, add P
        let is_neg = _mm256_cmpgt_epi64(zero, sum2);
        let res = _mm256_add_epi64(sum2, _mm256_and_si256(is_neg, p_val));

        _mm256_and_si256(res, mask_low)
    };

    let res_even = reduce_lane(even_64);
    let res_odd  = reduce_lane(odd_64);

    let odd_shifted = _mm256_slli_epi64::<32>(res_odd);
    _mm256_or_si256(res_even, odd_shifted)
}

// barrett_fake
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_7681(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(9));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(7681));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_10753(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(6));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(10753));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_11777(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(6));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(11777));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_12289(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(5));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(12289));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_13313(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(5));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(13313));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_15361(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(4));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(15361));
    let f = _mm256_sub_epi16(x, e);
    f
}

// """ Sage
// def to_int16(x):
//     u = Integer(x) & 0xFFFF 
//     if u >= 0x8000:
//         return u - 0x10000
//     return u

// f = [769, 3329, 7681, 7937, 9473, 10753]
// for i in f:
//     for j in range(10000):
//         a = ZZ.random_element(-i, i)
//         r = ZZ.random_element(-(2^15), 2^15-1)
//         m = ZZ.random_element(-(2^15), 2^15-1)
//         m1 = round((m*(2^16))/(i))
//         t = (r*m1) >> 16
//         d = to_int16(r*m)
//         dp = to_int16(a+d)
//         ds = to_int16(a-d)
//         ansp = to_int16(dp - to_int16(t*i))
//         anss = to_int16(ds + to_int16(t*i))
//         if ansp%i != (a + r*m)%i or ansp - (a + (r*m)%i) < -3*i or ansp - (a + (r*m)%i) > 3*i:
//             print(ansp, a + (r*m)%i, "||", a, r, m, ansp - (a + (r*m)%i))
//         if anss%i != (a - r*m)%i or anss - (a - (r*m)%i) < -3*i or anss - (a - (r*m)%i) > 3*i:
//             print(anss, a - (r*m)%i, "||", a, r, m, ansp - (a + (r*m)%i))
// """

// barrett_butterfly
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_7681(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_mullo_epi16(t, _mm256_set1_epi16(7681));
    let ansp = _mm256_sub_epi16(dp, e);
    let anss = _mm256_add_epi16(ds, e);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_10753(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_mullo_epi16(t, _mm256_set1_epi16(10753));
    let ansp = _mm256_sub_epi16(dp, e);
    let anss = _mm256_add_epi16(ds, e);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_11777(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_mullo_epi16(t, _mm256_set1_epi16(11777));
    let ansp = _mm256_sub_epi16(dp, e);
    let anss = _mm256_add_epi16(ds, e);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_12289(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_mullo_epi16(t, _mm256_set1_epi16(12289));
    let ansp = _mm256_sub_epi16(dp, e);
    let anss = _mm256_add_epi16(ds, e);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_13313(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_mullo_epi16(t, _mm256_set1_epi16(13313));
    let ansp = _mm256_sub_epi16(dp, e);
    let anss = _mm256_add_epi16(ds, e);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_15361(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_mullo_epi16(t, _mm256_set1_epi16(15361));
    let ansp = _mm256_sub_epi16(dp, e);
    let anss = _mm256_add_epi16(ds, e);
    [ansp, anss]
}


// barrett_fake_32 
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_7681_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(5585133);
    let m2 = _mm256_set1_epi32(7681);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(5585133));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(7681));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(5585133));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(7681));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
    ans
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_10753_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(1290167);
    let m2 = _mm256_set1_epi32(10753);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(1290167));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(10753));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(1290167));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(10753));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
    ans
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_11777_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(559168);
    let m2 = _mm256_set1_epi32(11777);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(559168));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(11777));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(559168));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(11777));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
    ans
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_12289_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(541132);
    let m2 = _mm256_set1_epi32(12289);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(541132));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(12289));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(541132));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(12289));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
    ans
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_13313_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(453390);
    let m2 = _mm256_set1_epi32(13313);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(453390));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(13313));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(453390));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(13313));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
    ans
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_15361_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(399420);
    let m2 = _mm256_set1_epi32(15361);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(399420));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(15361));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(399420));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(15361));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
    ans
}

// barrett_mul
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_7681(b: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(7681));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_10753(b: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(10753));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_11777(b: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(11777));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_12289(b: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(12289));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_13313(b: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(13313));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_15361(b: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(15361));
    let g = _mm256_sub_epi16(d, f);
    g
}

// barrett_ibutterfly
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_7681(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let ans0 = barrett_fake_7681(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul_7681(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_10753(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let ans0 = barrett_fake_10753(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul_10753(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_11777(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let ans0 = barrett_fake_11777(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul_11777(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_12289(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let ans0 = barrett_fake_12289(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul_12289(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_13313(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let ans0 = barrett_fake_13313(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul_13313(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_15361(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let ans0 = barrett_fake_15361(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul_15361(tmp, c, cr);
    [ans0, ans1]
}

// montproduct

#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_7681(x: __m256i, y: __m256i) -> __m256i{
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16(-7679));
    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16(7681));
    let f = _mm256_sub_epi16(hi, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_10753(x: __m256i, y: __m256i) -> __m256i{
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16(-10751));
    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16(10753));
    let f = _mm256_sub_epi16(hi, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_11777(x: __m256i, y: __m256i) -> __m256i{
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16(-11775));
    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16(11777));
    let f = _mm256_sub_epi16(hi, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_12289(x: __m256i, y: __m256i) -> __m256i{
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16(-12287));
    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16(12289));
    let f = _mm256_sub_epi16(hi, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_13313(x: __m256i, y: __m256i) -> __m256i{
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16(-13311));
    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16(13313));
    let f = _mm256_sub_epi16(hi, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_15361(x: __m256i, y: __m256i) -> __m256i{
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16(-15359));
    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16(15361));
    let f = _mm256_sub_epi16(hi, e);
    f
}


#[cfg(test)]
mod tests {
    use super::*; 
    use core::arch::x86_64::*;

    unsafe fn dump_m256i(v: __m256i) -> [i16; 16] {
        let mut arr = [0i16; 16];
        
        _mm256_storeu_si256(arr.as_mut_ptr() as *mut __m256i, v);
        arr
    }

    macro_rules! test_barrett {
        ($func_name:ident, $q:expr) => {
            #[test]
            fn $func_name() {
                if !is_x86_feature_detected!("avx2") {
                    return;
                }

                unsafe {
                    let inputs: [i16; 16] = [
                        0, 1, $q - 1, $q, 
                        $q + 1, 2 * $q, 2 * $q + 1,
                        -1, -10, -$q, 
                        32767, -32768, 
                        1024, 2048, 4096, 8192
                    ];

                    let input_vec = _mm256_loadu_si256(inputs.as_ptr() as *const __m256i);
                    let output_vec = super::$func_name(input_vec);
                    let outputs = dump_m256i(output_vec);

                    println!("\nTesting {} with q = {}", stringify!($func_name), $q);
                    println!("| Index |   Input |  Output | Diff (In - Out) | Status |");
                    println!("|-------|---------|---------|-----------------|--------|");

                    for i in 0..16 {
                        let x = inputs[i];
                        let y = outputs[i];
                        let diff = (x as i32) - (y as i32);
                        let is_valid = diff % ($q as i32) == 0;

                        println!("| {:5} | {:7} | {:7} | {:15} | {:^6} |", 
                            i, x, y, diff, if is_valid { "Ok" } else { "Failed" });

                        assert!(is_valid, 
                            "Reduction failed for {}: x={}, y={}, diff={}", 
                            stringify!($func_name), x, y, diff
                        );
                    }
                }
            }
        };
    }

    test_barrett!(barrett_fake_7681, 7681);
    test_barrett!(barrett_fake_10753, 10753);
    test_barrett!(barrett_fake_11777, 11777);
    test_barrett!(barrett_fake_12289, 12289);
    test_barrett!(barrett_fake_13313, 13313);
    test_barrett!(barrett_fake_15361, 15361);

    macro_rules! test_butterfly {
        ($func_name:ident, $q:expr) => {
            #[test]
            fn $func_name() {
                if !is_x86_feature_detected!("avx2") { return; }

                unsafe {
                    let q = $q as i32;

                    let c_test_cases = [1, q/2, -q/2, 2, 42];

                    let a_vals: [i16; 16] = [
                        0, 1, -1, (q-1) as i16, 
                        1000, -1000, 16384, -16384, 
                        (q/2) as i16, q as i16, -q as i16, 5,
                        0, 123, 456, 789
                    ];

                    let b_vals: [i16; 16] = [
                        0, 1, -1, (-1) as i16,
                        2, -2, 500, -500,
                        (q/4) as i16, 16384, -16384, 10,
                        q as i16, -q as i16, 1, -1
                    ];

                    for &c_val in &c_test_cases {
                        let c_i16 = c_val as i16;
                        let cr_val = ((c_i16 as f64 * 65536.0) / q as f64).round() as i16;

                        let a_vec = _mm256_loadu_si256(a_vals.as_ptr() as *const __m256i);
                        let b_vec = _mm256_loadu_si256(b_vals.as_ptr() as *const __m256i);
                        let c_vec = _mm256_set1_epi16(c_i16);
                        let cr_vec = _mm256_set1_epi16(cr_val);

                        let [ansp_vec, anss_vec] = super::$func_name(a_vec, b_vec, c_vec, cr_vec);
                        
                        let ansp_res = dump_m256i(ansp_vec);
                        let anss_res = dump_m256i(anss_vec);

                        for i in 0..16 {
                            let a = a_vals[i] as i32;
                            let b = b_vals[i] as i32;
                            let c = c_i16 as i32;

                            let expected_p = (a + b * c).rem_euclid(q);
                            let expected_s = (a - b * c).rem_euclid(q);

                            let res_p = (ansp_res[i] as i32).rem_euclid(q);
                            assert_eq!(res_p, expected_p, 
                                "\nFAILED P: q={}, c={}\nInput: a={}, b={}\nGot: {}, Expected: {}", 
                                q, c, a, b, ansp_res[i], expected_p);

                            let res_s = (anss_res[i] as i32).rem_euclid(q);
                            assert_eq!(res_s, expected_s, 
                                "\nFAILED S: q={}, c={}\nInput: a={}, b={}\nGot: {}, Expected: {}", 
                                q, c, a, b, anss_res[i], expected_s);
                        }
                    }
                    println!("SUCCESS: {} passed for all c cases.", stringify!($func_name));
                }
            }
        };
    }

    test_butterfly!(barrett_butterfly_7681, 7681);
    test_butterfly!(barrett_butterfly_10753, 10753);
    test_butterfly!(barrett_butterfly_11777, 11777);
    test_butterfly!(barrett_butterfly_12289, 12289);
    test_butterfly!(barrett_butterfly_13313, 13313);
    test_butterfly!(barrett_butterfly_15361, 15361);
    
    macro_rules! test_barrett_mul {
        ($func_name:ident, $q:expr) => {
            #[test]
            fn $func_name() {
                if !is_x86_feature_detected!("avx2") { return; }

                unsafe {
                    let q = $q as i32;

                    // 1. a within q/2 ~ -q/2
                    let test_a_factors = [1i16, 42, (q/2) as i16, (-q/2) as i16];
                    let b_vals: [i16; 16] = [
                        0, 1, -1, q as i16,
                        (q+1) as i16, (q-1) as i16, 1000, -1000,
                        (32767 / 7681 * 7681) as i16,
                        5, 10, 20, 50, 100, 200, 500
                    ];

                    for &a_val in &test_a_factors {
                        let a_i32 = a_val as i32;
                        let ar_overq_val = ((a_i32 as f64 * 65536.0) / q as f64).round() as i16;

                        let b_vec = _mm256_loadu_si256(b_vals.as_ptr() as *const __m256i);
                        let a_vec = _mm256_set1_epi16(a_val);
                        let ar_overq_vec = _mm256_set1_epi16(ar_overq_val);

                        let res_vec = super::$func_name(b_vec, a_vec, ar_overq_vec);
                        let results = dump_m256i(res_vec);

                        for i in 0..16 {
                            let b = b_vals[i] as i32;
                            let a = a_val as i32;
                            
                            let expected = (a * b).rem_euclid(q);

                            let actual = results[i] as i32;
                            let diff = (actual - expected).abs();

                            assert!(diff % q == 0, 
                                "\nFAILED: {} with q={}, a={}, ar={}\nInput b: {}\nGot: {}, Expected (mod q): {}\nDiff is not a multiple of q!", 
                                stringify!($func_name), q, a, ar_overq_val, b, actual, expected);
                            
                            assert!(actual > -2 * q && actual < 2 * q,
                                "Output range error: actual {}", actual);
                        }
                    }
                    println!("SUCCESS: {} passed.", stringify!($func_name));
                }
            }
        };
    }

    test_barrett_mul!(barrett_mul_7681, 7681);
    test_barrett_mul!(barrett_mul_10753, 10753);
    test_barrett_mul!(barrett_mul_11777, 11777);
    test_barrett_mul!(barrett_mul_12289, 12289);
    test_barrett_mul!(barrett_mul_13313, 13313);
    test_barrett_mul!(barrett_mul_15361, 15361);
}