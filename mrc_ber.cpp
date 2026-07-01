#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <filesystem>

unsigned long long nCr_under(int n_under, int r_under) {
    if (r_under < 0 || r_under > n_under) return 0;
    if (r_under == 0 || r_under == n_under) return 1;
    if (r_under > n_under / 2) r_under = n_under - r_under;
    
    unsigned long long res_under = 1;
    for (int k_under = 1; k_under <= r_under; ++k_under) {
        res_under = res_under * (n_under - k_under + 1) / k_under;
    }
    return res_under;
}

// BPSKの厳密式（等価SNR x_under を引数に取る）
double F_under(double x_under, int L_under) {
    double p_under = 0.5 * (1.0 - std::sqrt(x_under / (1.0 + x_under)));
    double Pb_under = 0.0;
    for (int k_under = 0; k_under < L_under; ++k_under) {
        double term_under = (double)nCr_under(L_under - 1 + k_under, k_under) * std::pow(1.0 - p_under, k_under);
        Pb_under += term_under;
    }
    Pb_under *= std::pow(p_under, L_under);
    return Pb_under;
}

double calculate_mrc_ber_under(double EbN0_lin_under, int L_under, int M_under) {
    if (M_under == 2) {
        return F_under(EbN0_lin_under, L_under);
    }
    
    // Square M-QAM (M = 4, 16, 64, 256...) の厳密BER式 (Gray coding)
    int sqM_under = std::round(std::sqrt(M_under));
    int log2sqM_under = std::round(std::log2(sqM_under));
    double log2M_under = std::log2(M_under);
    
    double sum_Pb_under = 0.0;
    
    for (int k_under = 1; k_under <= log2sqM_under; ++k_under) {
        int i_max_under = std::round((1.0 - std::pow(2.0, -k_under)) * sqM_under) - 1;
        for (int i_under = 0; i_under <= i_max_under; ++i_under) {
            
            double w_term1_under = std::floor( (i_under * std::pow(2.0, k_under - 1)) / sqM_under );
            int sign_under = ((int)w_term1_under % 2 == 1) ? -1 : 1;
            
            double w_term2_under = std::pow(2.0, k_under - 1) - std::floor( (i_under * std::pow(2.0, k_under - 1)) / sqM_under + 0.5 );
            double w_under = sign_under * w_term2_under;
            
            // 等価SNR x_under
            double x_under = (1.5 * std::pow(2 * i_under + 1, 2) * log2M_under / (M_under - 1.0)) * EbN0_lin_under;
            
            sum_Pb_under += w_under * F_under(x_under, L_under);
        }
    }
    
    double Pb_under = (2.0 / (sqM_under * log2sqM_under)) * sum_Pb_under;
    return Pb_under;
}

void main_under() {
    std::string out_dir_under = "C:\\Users\\Ide Nanako\\Desktop\\result2";
    std::filesystem::create_directories(out_dir_under);
    
    std::vector<int> M_values_under = {2, 4, 16, 64, 256};
    
    for (int L_under = 1; L_under <= 4; ++L_under) {
        std::string out_file_under = out_dir_under + "\\mrc_ber_results_L" + std::to_string(L_under) + ".csv";
        std::ofstream f_under(out_file_under);
        
        if (f_under.is_open()) {
            f_under << "EbN0_dB";
            for (int M_under : M_values_under) {
                if (M_under == 2) f_under << ",BPSK_BER";
                else f_under << "," << M_under << "QAM_BER";
            }
            f_under << "\n";
            
            for (int EbN0_dB_under = 0; EbN0_dB_under <= 30; EbN0_dB_under += 2) {
                double EbN0_lin_under = std::pow(10.0, EbN0_dB_under / 10.0);
                f_under << EbN0_dB_under;
                
                for (int M_under : M_values_under) {
                    double ber_under = calculate_mrc_ber_under(EbN0_lin_under, L_under, M_under);
                    f_under << "," << ber_under;
                }
                f_under << "\n";
            }
            f_under.close();
            std::cout << "L=" << L_under << " calculated and saved to: " << out_file_under << std::endl;
        } else {
            std::cerr << "Failed to open " << out_file_under << std::endl;
        }
    }
}

int main() {
    main_under();
    return 0;
}
