pub const TABLE_SIZE: usize = 4096;
pub const MAX_DIFF: f32 = 12.0;

// Lazy precomputed lookup table for ln(1 + exp(-delta))
pub struct LogSumTable {
    table: [f32; TABLE_SIZE],
    scale: f32,
}

impl LogSumTable {
    pub fn new() -> Self {
        let mut table = [0.0f32; TABLE_SIZE];
        let step = MAX_DIFF / (TABLE_SIZE as f32);
        for i in 0..TABLE_SIZE {
            let delta = (i as f32) * step;
            table[i] = (1.0f32 + (-delta).exp()).ln();
        }
        Self {
            table,
            scale: (TABLE_SIZE as f32) / MAX_DIFF,
        }
    }

    #[inline(always)]
    pub fn log_sum_exp(&self, a: f32, b: f32) -> f32 {
        if a.is_infinite() && a < 0.0 { return b; }
        if b.is_infinite() && b < 0.0 { return a; }

        let diff = (a - b).abs();
        let max_val = if a > b { a } else { b };

        if diff >= MAX_DIFF {
            max_val
        } else {
            let idx = (diff * self.scale) as usize;
            max_val + self.table[idx.min(TABLE_SIZE - 1)]
        }
    }
}

// Global static table
static GLOBAL_TABLE: std::sync::OnceLock<LogSumTable> = std::sync::OnceLock::new();

#[inline(always)]
pub fn log_sum_exp(a: f32, b: f32) -> f32 {
    GLOBAL_TABLE.get_or_init(LogSumTable::new).log_sum_exp(a, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_log_sum_exp() {
        let a = 2.0f32;
        let b = 3.0f32;
        let approx = log_sum_exp(a, b);
        let exact = (a.exp() + b.exp()).ln();
        assert!((approx - exact).abs() < 1e-3);
    }
}
