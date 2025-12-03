

pub fn optimize<F>(init_coeffs: &[f64], mut objective_fn: F) -> Vec<f64>
where
    F: FnMut(&[f64]) -> f64,
{
    let dim = init_coeffs.len();
    let alpha = 1.0;
    let gamma = 2.0;
    // Standard parameters
    let rho = 0.5;
    let sigma = 0.5;
    let max_iter = 1000; // Increased iteration count
    let tol = 1e-6;

    let mut simplex = Vec::with_capacity(dim + 1);
    simplex.push(init_coeffs.to_vec());

    for i in 0..dim {
        let mut p = init_coeffs.to_vec();
        // Better step size generation: proportional to value or fixed small step
        let step = if p[i].abs() > 1e-4 { p[i] * 0.05 } else { 0.01 };
        p[i] += step;
        simplex.push(p);
    }

    let mut scores: Vec<f64> = simplex.iter().map(|p| objective_fn(p)).collect();

    for _ in 0..max_iter {
        let mut indices: Vec<usize> = (0..dim + 1).collect();
        indices.sort_by(|&a, &b| scores[a].partial_cmp(&scores[b]).unwrap_or(std::cmp::Ordering::Equal));

        let best_idx = indices[0];
        let worst_idx = indices[dim];
        let second_worst_idx = indices[dim - 1];

        let best_score = scores[best_idx];
        let worst_score = scores[worst_idx];

        if (worst_score - best_score).abs() < tol {
            break;
        }

        let mut centroid = vec![0.0; dim];
        for i in 0..dim { // exclude worst
            let idx = indices[i];
            for d in 0..dim {
                centroid[d] += simplex[idx][d];
            }
        }
        for d in 0..dim {
            centroid[d] /= dim as f64;
        }

        // Reflection
        let diff_c_w = sub(&centroid, &simplex[worst_idx]);
        let xr = add(&centroid, &scale(&diff_c_w, alpha));
        let score_r = objective_fn(&xr);

        if score_r < scores[second_worst_idx] && score_r >= best_score {
            simplex[worst_idx] = xr;
            scores[worst_idx] = score_r;
            continue;
        }

        // Expansion
        if score_r < best_score {
            let diff_xr_c = sub(&xr, &centroid);
            let xe = add(&centroid, &scale(&diff_xr_c, gamma));
            let score_e = objective_fn(&xe);
            if score_e < score_r {
                simplex[worst_idx] = xe;
                scores[worst_idx] = score_e;
            } else {
                simplex[worst_idx] = xr;
                scores[worst_idx] = score_r;
            }
            continue;
        }

        // Contraction
        let diff_w_c = sub(&simplex[worst_idx], &centroid);
        let xc = add(&centroid, &scale(&diff_w_c, rho));
        let score_c = objective_fn(&xc);
        
        if score_c < scores[worst_idx] {
            simplex[worst_idx] = xc;
            scores[worst_idx] = score_c;
            continue;
        }

        // Shrink
        for i in 1..=dim {
            let idx = indices[i];
            let diff = sub(&simplex[idx], &simplex[best_idx]);
            simplex[idx] = add(&simplex[best_idx], &scale(&diff, sigma));
            scores[idx] = objective_fn(&simplex[idx]);
        }
    }

    let best_idx = indices_of_min(&scores);
    simplex[best_idx].clone()
}

// Helper math functions
fn add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b).map(|(x, y)| x + y).collect()
}

fn sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b).map(|(x, y)| x - y).collect()
}

fn scale(a: &[f64], s: f64) -> Vec<f64> {
    a.iter().map(|x| x * s).collect()
}

fn indices_of_min(v: &[f64]) -> usize {
    v.iter().enumerate().min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)).map(|(i, _)| i).unwrap()
}
