// https://www.codingame.com/ide/puzzle/running-up-that-hill

use std::io;

macro_rules! parse_input {
    ($x:expr, $t:ident) => ($x.trim().parse::<$t>().unwrap())
}

// fast optimized hill cipher solver for modulus 45 using inv mod5 and mod9 and crt
// gauss jordan used for inversion mod 5 and mod 9 for speed

fn build_charset() -> (Vec<char>, [i64; 256]) {
    // charset and lookup table
    let s = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:";
    let mut charset: Vec<char> = s.chars().collect();
    let mut lut: [i64; 256] = [-1i64; 256];
    for (i, &ch) in charset.iter().enumerate() {
        let idx = ch as usize;
        lut[idx] = i as i64;
    }
    (charset, lut)
}

fn mod_norm_i64(x: i64, m: i64) -> i64 {
    let mut r = x % m;
    if r < 0 { r += m; }
    r
}

fn egcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a.abs(), if a >= 0 { 1 } else { -1 }, 0)
    } else {
        let (g, x1, y1) = egcd(b, a % b);
        (g, y1, x1 - (a / b) * y1)
    }
}

fn inv_mod_any(a: i64, modu: i64) -> Option<i64> {
    let (g, x, _) = egcd(a, modu);
    if g != 1 {
        None
    } else {
        Some(mod_norm_i64(x, modu))
    }
}

// gauss jordan inverse for modulus m works for m=5 or m=9 here
// returns none if matrix not invertible mod m
fn mat_inv_mod(mut a: Vec<Vec<i64>>, n: usize, modu: i64) -> Option<Vec<Vec<i64>>> {
    // augment with identity matrix on the right
    let mut aug = vec![vec![0i64; 2 * n]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = mod_norm_i64(a[i][j], modu);
        }
        for j in 0..n {
            aug[i][n + j] = if i == j { 1 } else { 0 };
        }
    }

    for col in 0..n {
        // find pivot row with element invertible modulo modu
        let mut pivot_row = None;
        for r in col..n {
            let val = aug[r][col];
            if egcd(val, modu).0 == 1 {
                pivot_row = Some(r);
                break;
            }
        }
        if pivot_row.is_none() {
            return None;
        }
        let pr = pivot_row.unwrap();
        if pr != col {
            aug.swap(pr, col);
        }
        let pivot = aug[col][col];
        let inv = inv_mod_any(pivot, modu).unwrap();
        // normalize pivot row
        for j in col..2 * n {
            aug[col][j] = mod_norm_i64(aug[col][j] * inv, modu);
        }
        // eliminate other rows
        for r in 0..n {
            if r == col { continue; }
            let factor = aug[r][col];
            if factor != 0 {
                for j in col..2 * n {
                    let v = aug[r][j] - factor * aug[col][j];
                    aug[r][j] = mod_norm_i64(v, modu);
                }
            }
        }
    }
    // extract inverse from right half
    let mut inv = vec![vec![0i64; n]; n];
    for i in 0..n {
        for j in 0..n {
            inv[i][j] = mod_norm_i64(aug[i][n + j], modu);
        }
    }
    Some(inv)
}

// combine inverse matrices modulo 5 and 9 into modulo 45 using crt
// solve x ≡ a9 (mod9) and x ≡ a5 (mod5)
// formula uses t = 4 * (a5 - a9) mod5 then x = a9 + 9 * t
fn combine_inv_mod45(inv5: &Vec<Vec<i64>>, inv9: &Vec<Vec<i64>>, n: usize) -> Vec<Vec<i64>> {
    let mut res = vec![vec![0i64; n]; n];
    for i in 0..n {
        for j in 0..n {
            let a5 = mod_norm_i64(inv5[i][j], 5);
            let a9 = mod_norm_i64(inv9[i][j], 9);
            let rhs = mod_norm_i64(a5 - a9, 5);
            let t = mod_norm_i64(4 * rhs, 5);
            let x = a9 + 9 * t;
            res[i][j] = mod_norm_i64(x, 45);
        }
    }
    res
}

// multiply matrices modulo modu
fn mat_mul_mod(a: &Vec<Vec<i64>>, b: &Vec<Vec<i64>>, modu: i64) -> Vec<Vec<i64>> {
    let r = a.len();
    let m = a[0].len();
    let c = b[0].len();
    let mut res = vec![vec![0i64; c]; r];
    for i in 0..r {
        let ai = &a[i];
        let ri = &mut res[i];
        for k in 0..m {
            let av = ai[k];
            if av == 0 { continue; }
            let bv = &b[k];
            for j in 0..c {
                let tmp = ri[j] + av * bv[j];
                ri[j] = if tmp >= 0 {
                    tmp % modu
                } else {
                    ((tmp % modu) + modu) % modu
                };
            }
        }
    }
    for i in 0..r {
        for j in 0..c {
            res[i][j] = mod_norm_i64(res[i][j], modu);
        }
    }
    res
}

// generate combinations of m indices choose k
fn combinations_indices(m: usize, k: usize) -> Vec<Vec<usize>> {
    let mut res = Vec::new();
    let mut comb = Vec::with_capacity(k);
    fn backtrack(start: usize, m: usize, k: usize, comb: &mut Vec<usize>, res: &mut Vec<Vec<usize>>) {
        if comb.len() == k {
            res.push(comb.clone());
            return;
        }
        for i in start..m {
            comb.push(i);
            backtrack(i + 1, m, k, comb, res);
            comb.pop();
        }
    }
    backtrack(0, m, k, &mut comb, &mut res);
    res
}

// fast text to numbers using lookup table
fn text_to_nums_fast(s: &str, lut: &[i64; 256]) -> Vec<i64> {
    let mut v = Vec::with_capacity(s.len());
    for ch in s.bytes() {
        let idx = lut[ch as usize];
        v.push(idx);
    }
    v
}

fn nums_to_text(v: &Vec<i64>, charset: &Vec<char>) -> String {
    let mut out = String::with_capacity(v.len());
    for &x in v.iter() {
        let idx = (x.rem_euclid(45)) as usize;
        out.push(charset[idx]);
    }
    out
}

fn divisors_of(n: usize) -> Vec<usize> {
    let mut d = Vec::new();
    for i in 2..=n {
        if n % i == 0 {
            d.push(i);
        }
    }
    d
}

fn main() {
    let mut input_line = String::new();
    io::stdin().read_line(&mut input_line).unwrap();
    let cipher = input_line.trim_end_matches('\n').to_string();
    input_line.clear();
    io::stdin().read_line(&mut input_line).unwrap();
    let clear = input_line.trim_end_matches('\n').to_string();
    input_line.clear();
    io::stdin().read_line(&mut input_line).unwrap();
    let cipher_to_decipher = input_line.trim_end_matches('\n').to_string();
    input_line.clear();
    io::stdin().read_line(&mut input_line).unwrap();
    let clear_to_cipher = input_line.trim_end_matches('\n').to_string();

    let (charset, lut) = build_charset();

    let modu = 45i64;

    let clear_nums = text_to_nums_fast(&clear, &lut);
    let cipher_nums = text_to_nums_fast(&cipher, &lut);

    let L = clear_nums.len();
    let divisors = divisors_of(L);

    let mut found = false;
    let mut best_n = 0usize;
    let mut best_A: Vec<Vec<i64>> = Vec::new();

    'outer: for &n in divisors.iter() {
        let mblocks = L / n;
        let mut X = vec![vec![0i64; mblocks]; n];
        let mut Y = vec![vec![0i64; mblocks]; n];
        for b in 0..mblocks {
            for i in 0..n {
                X[i][b] = clear_nums[b * n + i];
                Y[i][b] = cipher_nums[b * n + i];
            }
        }

        let combs = combinations_indices(mblocks, n);
        for comb in combs {
            let mut X0 = vec![vec![0i64; n]; n];
            let mut Y0 = vec![vec![0i64; n]; n];
            for col in 0..n {
                let bidx = comb[col];
                for row in 0..n {
                    X0[row][col] = X[row][bidx];
                    Y0[row][col] = Y[row][bidx];
                }
            }
            // try invert X0 mod5 and mod9 and combine
            if let (Some(inv5), Some(inv9)) = (mat_inv_mod(X0.clone(), n, 5), mat_inv_mod(X0.clone(), n, 9)) {
                let inv45 = combine_inv_mod45(&inv5, &inv9, n);
                let A = mat_mul_mod(&Y0, &inv45, modu);
                let prod = mat_mul_mod(&A, &X, modu);
                let mut ok = true;
                for i in 0..n {
                    for j in 0..mblocks {
                        if mod_norm_i64(prod[i][j], modu) != mod_norm_i64(Y[i][j], modu) {
                            ok = false;
                            break;
                        }
                    }
                    if !ok { break; }
                }
                if ok {
                    found = true;
                    best_n = n;
                    best_A = A;
                    break 'outer;
                }
            }
        }
    }

    if !found {
        eprintln!("unable to find an invertible block set");
        // print empty outputs to match required format
        println!();
        println!();
        return;
    }

    let n = best_n;
    let invA5 = mat_inv_mod(best_A.clone(), n, 5).expect("a invertible mod5");
    let invA9 = mat_inv_mod(best_A.clone(), n, 9).expect("a invertible mod9");
    let invA = combine_inv_mod45(&invA5, &invA9, n);

    // decipher line 3
    let ctd_nums = text_to_nums_fast(&cipher_to_decipher, &lut);
    let mblocks_ctd = ctd_nums.len() / n;
    let mut Cmat = vec![vec![0i64; mblocks_ctd]; n];
    for b in 0..mblocks_ctd {
        for i in 0..n {
            Cmat[i][b] = ctd_nums[b * n + i];
        }
    }
    let plain_mat = mat_mul_mod(&invA, &Cmat, modu);
    let mut plain_nums: Vec<i64> = Vec::with_capacity(mblocks_ctd * n);
    for b in 0..mblocks_ctd {
        for i in 0..n {
            plain_nums.push(mod_norm_i64(plain_mat[i][b], modu));
        }
    }
    let plaintext = nums_to_text(&plain_nums, &charset);

    // cipher line 4
    let ctc_nums = text_to_nums_fast(&clear_to_cipher, &lut);
    let mblocks_ctc = ctc_nums.len() / n;
    let mut Pmat = vec![vec![0i64; mblocks_ctc]; n];
    for b in 0..mblocks_ctc {
        for i in 0..n {
            Pmat[i][b] = ctc_nums[b * n + i];
        }
    }
    let ciphered_mat = mat_mul_mod(&best_A, &Pmat, modu);
    let mut ciphered_nums: Vec<i64> = Vec::with_capacity(mblocks_ctc * n);
    for b in 0..mblocks_ctc {
        for i in 0..n {
            ciphered_nums.push(mod_norm_i64(ciphered_mat[i][b], modu));
        }
    }
    let ciphered = nums_to_text(&ciphered_nums, &charset);

    println!("{}", plaintext);
    println!("{}", ciphered);
}
