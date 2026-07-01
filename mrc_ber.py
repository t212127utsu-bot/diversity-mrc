import os
import math
import csv

def nCr_under(n_under, r_under):
    if r_under < 0 or r_under > n_under:
        return 0
    try:
        return math.comb(n_under, r_under)
    except AttributeError:
        f = math.factorial
        return f(n_under) // f(r_under) // f(n_under - r_under)

def calculate_mrc_ber_under(EbN0_lin_under, L_under):
    p_under = 0.5 * (1.0 - math.sqrt(EbN0_lin_under / (1.0 + EbN0_lin_under)))
    
    Pb_under = 0.0
    for k_under in range(L_under):
        term_under = nCr_under(L_under - 1 + k_under, k_under) * ( (1.0 - p_under) ** k_under )
        Pb_under += term_under
        
    Pb_under *= (p_under ** L_under)
    return Pb_under

def main_under():
    # Output directory
    out_dir_under = r"C:\Users\Ide Nanako\Desktop\result2"
    if not os.path.exists(out_dir_under):
        os.makedirs(out_dir_under)
        
    # Dictionary to store results for each L
    results_under = {1: [], 2: [], 3: [], 4: []}
    
    # Calculate BER for EbN0 from 0 to 30 dB (step 2)
    for EbN0_dB_under in range(0, 31, 2):
        EbN0_lin_under = 10.0 ** (EbN0_dB_under / 10.0)
        
        for L_under in range(1, 5):
            ber_under = calculate_mrc_ber_under(EbN0_lin_under, L_under)
            results_under[L_under].append([EbN0_dB_under, ber_under])
            
    # Write results to separate CSV files
    for L_under in range(1, 5):
        out_file_under = os.path.join(out_dir_under, f"l{L_under}_mrc_results.csv")
        with open(out_file_under, mode='w', newline='') as f_under:
            writer_under = csv.writer(f_under)
            header_under = ["EbN0_dB", f"L={L_under}_BER"]
            writer_under.writerow(header_under)
            
            for row_under in results_under[L_under]:
                writer_under.writerow(row_under)
                
        print(f"L={L_under}の計算完了。結果を保存しました: {out_file_under}")

if __name__ == "__main__":
    main_under()
