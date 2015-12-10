functions {
  int first_capture(int[] y_i) {
    for (k in 1:size(y_i))
      if (y_i[k])
        return k;
    return 0;
  }
  int last_capture(int[] y_i) {
    for (k_rev in 0:(size(y_i) - 1)) {
      int k;
      k <- size(y_i) - k_rev;
      if (y_i[k])
        return k;
    }
    return 0;
  }
  real state_space_log(int[] y, vector phi, vector p, vector gamma) {
    int ft;
    int lt;
    int S;
    vector[first_capture(y) * (size(y) - last_capture(y) + 1)] lp;
    int i;

    ft <- first_capture(y);
    lt <- last_capture(y);
    S <- size(y);
    i <- 1;

    for(t_first_alive in 1:ft) {
      for (t_last_alive in lt:S) {
        real sl;
        int z[S];

        for(l in 1:S) {
          z[l] <- 0;
        }
        for(a in t_first_alive:t_last_alive) {
          z[a] <- 1;
        }

        sl <- bernoulli_log(z[1], gamma[1]);
        for(j in 2:S) {
          sl <- sl + bernoulli_log(z[j], (z[j - 1] * phi[j - 1]) + 
              (gamma[j] * (1 - z[j - 1])));
        }
        for(k in 1:S) {
          sl <- sl + bernoulli_log(y[k], z[k] * p[k]);
        }

        lp[i] <- sl;
        i <- i + 1;
      }
    }
    return log_sum_exp(lp);
  }
}
data {
  int N;
  int T;
  int sight[N, T];
}
parameters {
  vector<lower=0,upper=1>[T - 1] phi;
  vector<lower=0,upper=1>[T] p;
  vector<lower=0,upper=1>[T] gamma;
}
model {
  for(n in 1:N) {
    sight[n] ~ state_space(phi, p, gamma);
  }
}
