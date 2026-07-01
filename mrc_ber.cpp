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

double calculate_mrc_ber_under(double EbN0_lin_under, int L_under) {
    double p_under = 0.5 * (1.0 - std::sqrt(EbN0_lin_under / (1.0 + EbN0_lin_under)));
    
    double Pb_under = 0.0;
    for (int k_under = 0; k_under < L_under; ++k_under) {
        double term_under = (double)nCr_under(L_under - 1 + k_under, k_under) * std::pow(1.0 - p_under, k_under);
        Pb_under += term_under;
    }
    
    Pb_under *= std::pow(p_under, L_under);
    return Pb_under;
}

void main_under() {
    std::string out_dir_under = "C:\\Users\\Ide Nanako\\Desktop\\result2";
    std::filesystem::create_directories(out_dir_under);
    
    std::vector<std::vector<std::pair<int, double>>> results_under(5);
    
    for (int EbN0_dB_under = 0; EbN0_dB_under <= 30; EbN0_dB_under += 2) {
        double EbN0_lin_under = std::pow(10.0, EbN0_dB_under / 10.0);
        
        for (int L_under = 1; L_under <= 4; ++L_under) {
            double ber_under = calculate_mrc_ber_under(EbN0_lin_under, L_under);
            results_under[L_under].push_back({EbN0_dB_under, ber_under});
        }
    }
    
    for (int L_under = 1; L_under <= 4; ++L_under) {
        std::string out_file_under = out_dir_under + "\\l" + std::to_string(L_under) + "_mrc_results.csv";
        std::ofstream f_under(out_file_under);
        if (f_under.is_open()) {
            f_under << "EbN0_dB,L=" << L_under << "_BER\n";
            for (const auto& row_under : results_under[L_under]) {
                f_under << row_under.first << "," << row_under.second << "\n";
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
