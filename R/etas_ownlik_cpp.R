# C++ version of the likelihood in etas_ownlik.R (about 50x faster), used for the bootstrap.
# Copyright (C) 2026  Davood Kazemi Lafmejani and Elham Shabani. GPL-3.0-or-later.
Rcpp::cppFunction('
double ll_cpp(NumericVector par, NumericVector t, NumericVector x, NumericVector y, NumericVector m,
              LogicalVector flag, NumericVector bk, double integ0, double T1, NumericVector xr, NumericVector yr,
              double r, double tb) {
  int N = t.size();
  double mu = par[0], A = par[1], c = par[2], al = par[3], p = par[4], D = par[5], ga = par[6];
  double ll = 0, comp_trig = 0;
  std::vector<double> kap(N), sg(N);
  for (int i = 0; i < N; ++i) { kap[i] = A * std::exp(al * m[i]); sg[i] = D * std::exp(ga * m[i]); }
  for (int j = 0; j < N; ++j) {
    if (flag[j]) {
      double s = mu * bk[j];
      for (int i = 0; i < j; ++i) {
        double dt = t[j] - t[i];
        if (dt <= 0) continue;
        double dx = x[j] - x[i], dy = y[j] - y[i], s2 = sg[i] * sg[i];
        s += kap[i] * (p - 1) / c * std::pow(1 + dt / c, -p) * std::exp(-(dx * dx + dy * dy) / (2 * s2)) / (2 * M_PI * s2);
      }
      double rj = (t[j] > tb) ? r : 1.0;
      ll += std::log(std::max(rj * s, 1e-25));
    }
    double sj = (R::pnorm((xr[1] - x[j]) / sg[j], 0, 1, 1, 0) - R::pnorm((xr[0] - x[j]) / sg[j], 0, 1, 1, 0)) *
                (R::pnorm((yr[1] - y[j]) / sg[j], 0, 1, 1, 0) - R::pnorm((yr[0] - y[j]) / sg[j], 0, 1, 1, 0));
    double Gl = 1 - std::pow(1 + (T1 - t[j]) / c, -(p - 1));
    double Gt = (tb < T1 && tb > t[j]) ? 1 - std::pow(1 + (tb - t[j]) / c, -(p - 1)) : ((tb <= t[j]) ? 0.0 : Gl);
    if (!(tb < T1)) comp_trig += kap[j] * sj * Gl;
    else comp_trig += kap[j] * sj * (Gt + r * (Gl - Gt));
  }
  double comp_bg = (tb < T1) ? mu * integ0 * (tb / T1 + r * (1 - tb / T1)) : mu * integ0;
  return ll - comp_bg - comp_trig;
}')
