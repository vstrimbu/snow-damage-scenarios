// --- SNOW DAMAGE FUNCTIONS FOR GAYA 2.0

// TODO: the birch coefficients are taken from pine
static vector<vector<double>> sb_coef = {
    {36.0e+06, 46.0e+06, 46.0e+06}, // 0. MOR (MPA)
    {0.90, 0.85, 0.85}, // 1. fknot
    {6.30e+09, 7.3e+09, 7.3e+09}, // 2. MOE (N*m^-2)
    {128.5, 142.9, 142.9}, // 3. Creg (N*m*kg^-1)
    {0.000068, 0.000086, 0.000086}, // 4. SWparam0 (m^3)
    {1.827421, 1.634205, 1.634205}, // 5. SWparam1 
    {0.981362, 1.079539, 1.079539}, // 6. SWparam2 
    {850, 850, 850}, // 7. stemdens (kg*m^-3)
    {0.4206, 1.3424, 1.3424}, // 8. cdepthparam0 (m)
    {0.4368, 0.3156, 0.3156}, // 9. cdepthparam1
    {0.5824, 0.7005, 0.7005}, // 10. cwidthparam0 (m)
    {0.11500, 0.11980, 0.11980}, // 11. cwidthparam1
    {2.50, 2.50, 2.50}, // 12. cdensity (kg*m^-3)
};

const double g = 9.81; // Earth's gravitational acceleration (m*s^-2)

/**
 * Predict tree crown depth
 * @param sp Species (1 = spruce, 2 = pine, 3 = birch) 
 * @param h Height (m) 
 * @return Crown depth (m) 
 */
double predict_crown_depth(double sp, double h) {
    return sb_coef[8][sp - 1] + sb_coef[9][sp - 1] * h;
}

/**
 * Predict tree crown width
 * @param sp Species (1 = spruce, 2 = pine, 3 = birch) 
 * @param d Diameter at breast height (cm) 
 * @return Crown width (m) 
 */
double predict_crown_width(double sp, double d) {
    return sb_coef[10][sp - 1] + sb_coef[11][sp - 1] * d;
}

/**
 * Predict tree crown weight
 * @param sp Species (1 = spruce, 2 = pine, 3 = birch)
 * @param h Height (m) 
 * @param d Diameter at breast height (cm) 
 * @return Crown weight (kg) 
 */
double predict_crown_weight(double sp, double h, double d) {
    return (M_PI * sb_coef[12][sp - 1] / 3) * predict_crown_depth(sp, h) * pow(predict_crown_width(sp, d) / 2, 2);
}

/**
 * Predict critical snow load 
 * @param sp Species (1 = spruce, 2 = pine, 3 = birch)
 * @param h Height (m) 
 * @param d Diameter at breast height (cm) 
 * @return Critical snow load per tree (kg/m2)
 * @reference Zubkov et al. 2023
 * @reference Locatelli et al. 2022
 */
double predict_critical_snow_load(double sp, double h, double d) {
    double cw = predict_crown_weight(sp, h, d);
    double a = predict_crown_depth(sp, h) / (2 * h);
    double result = 1 / (((64 * pow(h, 2) * g) / (M_PI * sb_coef[2][sp - 1] * pow(d / 100, 4)))*(5 * a - 5.96 * pow(a, 0.6) - 0.71 * pow(a, 2) + 1.67)) - cw;
    double crown_area = pow(0.5 * predict_crown_width(sp, d), 2) * M_PI;
    crown_area = crown_area <= 0 ? 1 : crown_area;
    return result > 0 ? result / crown_area : 0.01;
}

/**
 * Predict tree breaking probability
 * @param csl Critical snow load (kg)
 * @param asl Actual snow load (kg) 
 * @return Breaking probability (0..1)
 * @reference Zubkov et al. 2023 
 */
double predict_breaking_probability(double csl, double asl) {
    return 1 / (1 + exp((csl - asl) / 1000));
}

/**
 * Predict tree breaking probability
 * @param sp Species (1 = spruce, 2 = pine, 3 = birch)
 * @param h Height (m) 
 * @param d Diameter at breast height (cm)  
 * @param asl Actual snow load (kg) 
 * @return Breaking probability (0..1)
 * @reference Zubkov et al. 2023 
 */
double predict_breaking_probability(double sp, double h, double d, double asl) {
    double csl = predict_critical_snow_load(sp, h, d);
    return predict_breaking_probability(csl, asl);
}