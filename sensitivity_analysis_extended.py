# -*- coding: utf-8 -*-
"""
擴展敏感度分析腳本 (Extended Sensitivity Analysis)
===================================================
版本: v2.0
日期: 2025-12-17

包含: UVA逆境、花青素、碳依賴修復、LDMC 共 25 個參數
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from scipy.integrate import solve_ivp

from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment, ODE_SETTINGS

# 擴展的參數列表 - 包含花青素和碳依賴修復
PARAMS_TO_ANALYZE = {
    # === UVA 逆境機制 (13個) ===
    'E_50': {'base': 237.6, 'range': 0.30, 'cat': 'UVA_Stress', 'desc': '半效能量'},
    'E_scale': {'base': 39.6, 'range': 0.30, 'cat': 'UVA_Stress', 'desc': '能量尺度'},
    'stress_damage_coeff': {'base': 3.5e-6, 'range': 0.30, 'cat': 'UVA_Stress', 'desc': '損傷係數'},
    'stress_repair_coeff': {'base': 1.0e-5, 'range': 0.30, 'cat': 'UVA_Stress', 'desc': '修復係數'},
    'LAI_ref_vuln': {'base': 7.5, 'range': 0.30, 'cat': 'UVA_Stress', 'desc': 'LAI參考值'},
    'n_vuln': {'base': 7, 'range': 0.30, 'cat': 'UVA_Stress', 'desc': '脆弱性指數'},
    'k_intraday': {'base': 1.5, 'range': 0.30, 'cat': 'UVA_Stress', 'desc': '日內非線性'},
    'stress_nonlinear_coeff': {'base': 1.5, 'range': 0.30, 'cat': 'UVA_Stress', 'desc': '非線性係數'},
    'K_nonlinear': {'base': 3.0, 'range': 0.30, 'cat': 'UVA_Stress', 'desc': '非線性半飽和'},
    'circadian_disruption_factor': {'base': 2.0, 'range': 0.50, 'cat': 'UVA_Stress', 'desc': '節律損傷'},
    'stress_photosynthesis_inhibition': {'base': 0.70, 'range': 0.30, 'cat': 'UVA_Stress', 'desc': '光合抑制'},
    'K_stress': {'base': 5.0, 'range': 0.30, 'cat': 'UVA_Stress', 'desc': '抑制半飽和'},
    'par_conversion_factor': {'base': 1.0, 'range': 0.10, 'cat': 'UVA_Stress', 'desc': 'UVA-PAR轉換'},

    # === 花青素機制 (7個) ===
    'base_anth_rate_light': {'base': 4.11e-10, 'range': 0.30, 'cat': 'Anthocyanin', 'desc': '日間基礎合成'},
    'base_anth_rate_dark': {'base': 2.05e-10, 'range': 0.30, 'cat': 'Anthocyanin', 'desc': '夜間基礎合成'},
    'V_max_anth': {'base': 2.5e-10, 'range': 0.30, 'cat': 'Anthocyanin', 'desc': '最大誘導合成'},
    'K_m_anth': {'base': 30.0, 'range': 0.30, 'cat': 'Anthocyanin', 'desc': '半飽和Stress'},
    'stress_threshold_anth': {'base': 15.0, 'range': 0.30, 'cat': 'Anthocyanin', 'desc': '誘導閾值'},
    'anth_threshold_sharpness': {'base': 0.5, 'range': 0.30, 'cat': 'Anthocyanin', 'desc': '閾值銳度'},
    'k_deg': {'base': 2.5e-6, 'range': 0.30, 'cat': 'Anthocyanin', 'desc': '降解率'},

    # === 碳依賴修復 (3個) ===
    'base_repair_capacity': {'base': 0.5, 'range': 0.30, 'cat': 'Carbon_Repair', 'desc': '基礎修復'},
    'carbon_repair_bonus': {'base': 0.5, 'range': 0.30, 'cat': 'Carbon_Repair', 'desc': '碳池加成'},
    'K_carbon': {'base': 0.001, 'range': 0.30, 'cat': 'Carbon_Repair', 'desc': '碳池半飽和'},

    # === LDMC (2個) ===
    'ldmc_stress_sensitivity': {'base': 1.0, 'range': 0.30, 'cat': 'LDMC', 'desc': 'LDMC敏感度'},
    'K_ldmc': {'base': 50.0, 'range': 0.30, 'cat': 'LDMC', 'desc': 'LDMC半飽和'},
}

TEST_TREATMENTS = ['CK', 'L6D6', 'L6D6-N', 'H12D3']

def run_sim(p, treatment):
    """執行單次模擬，返回 FW, Stress, Anth_ppm"""
    env = get_env_for_treatment(treatment)
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    initial_state = [Xd_init, 0.0, LAI_init, Anth_init, 0.0]

    try:
        sol = solve_ivp(
            uva_sun_derivatives,
            (0, SIMULATION['days'] * 86400),
            initial_state,
            args=(p, env),
            method=ODE_SETTINGS['method'],
            max_step=ODE_SETTINGS['max_step']
        )
        if sol.success:
            stress = sol.y[4, -1]
            anth = sol.y[3, -1]
            ratio = calculate_dynamic_dw_fw_ratio(stress, p)
            dw = sol.y[0, -1] / ENV_BASE['plant_density'] * 1000
            fw = dw / ratio
            # 計算花青素 ppm
            fw_kg_m2 = fw * ENV_BASE['plant_density'] / 1000
            anth_ppm = anth * 1e6 / (fw_kg_m2 + 1e-9)
            return fw, stress, anth_ppm, True
    except Exception:
        pass
    return np.nan, np.nan, np.nan, False


def main():
    print('=' * 80)
    print('擴展敏感度分析 - 萵苣UVA模型 v5.7 (含花青素與碳修復)')
    print('=' * 80)
    print()
    print(f'分析參數數量: {len(PARAMS_TO_ANALYZE)}')
    print(f'測試處理組: {TEST_TREATMENTS}')
    print()

    # 按類別分組顯示
    categories = {}
    for pname, pinfo in PARAMS_TO_ANALYZE.items():
        cat = pinfo['cat']
        if cat not in categories:
            categories[cat] = []
        categories[cat].append(pname)

    for cat, params in categories.items():
        print(f'  {cat}: {len(params)} 個參數')
    print()

    # 執行分析
    all_results = {}
    total = len(PARAMS_TO_ANALYZE)

    for i, (pname, pinfo) in enumerate(PARAMS_TO_ANALYZE.items()):
        base = pinfo['base']
        var = pinfo['range']
        cat = pinfo['cat']

        print(f'分析 [{i+1}/{total}]: {pname}...', end=' ', flush=True)

        results = {'treatments': {}}
        for treatment in TEST_TREATMENTS:
            p_base = UVAParams()
            fw_b, st_b, anth_b, ok_b = run_sim(p_base, treatment)

            p_low = UVAParams()
            setattr(p_low, pname, base * (1 - var))
            fw_l, st_l, anth_l, ok_l = run_sim(p_low, treatment)

            p_high = UVAParams()
            setattr(p_high, pname, base * (1 + var))
            fw_h, st_h, anth_h, ok_h = run_sim(p_high, treatment)

            if ok_b and ok_l and ok_h:
                dx = 2 * var
                # FW 敏感度
                dfw = (fw_h - fw_l) / fw_b if fw_b > 0 else 0
                sens_fw = dfw / dx
                # Stress 敏感度
                dst = (st_h - st_l) / (st_b + 0.01)
                sens_st = dst / dx
                # Anth 敏感度
                danth = (anth_h - anth_l) / (anth_b + 0.01)
                sens_anth = danth / dx

                results['treatments'][treatment] = {
                    'FW_sens': sens_fw,
                    'Stress_sens': sens_st,
                    'Anth_sens': sens_anth,
                    'FW_base': fw_b,
                    'Anth_base': anth_b,
                    'Stress_base': st_b,
                }

        all_results[pname] = {'cat': cat, 'results': results}
        print('Done')

    # 打印結果
    print()
    print('=' * 100)
    print('敏感度分析結果 (彈性係數 S = dY/Y / dX/X)')
    print('=' * 100)

    for treatment in TEST_TREATMENTS:
        print()
        print(f'[{treatment}]')
        print('-' * 95)
        header = f'{"Parameter":<32} {"Category":<15} {"FW":>10} {"Stress":>10} {"Anth":>10}'
        print(header)
        print('-' * 95)

        sens_list = []
        for pname, data in all_results.items():
            cat = data['cat']
            if treatment in data['results']['treatments']:
                r = data['results']['treatments'][treatment]
                fw_s = r['FW_sens']
                st_s = r['Stress_sens']
                anth_s = r['Anth_sens']
                row = f'{pname:<32} {cat:<15} {fw_s:>10.3f} {st_s:>10.3f} {anth_s:>10.3f}'
                print(row)
                sens_list.append((pname, cat, abs(fw_s), abs(st_s), abs(anth_s)))

        print()
        # Top 3 for each metric
        top3_fw = sorted(sens_list, key=lambda x: x[2], reverse=True)[:3]
        top3_st = sorted(sens_list, key=lambda x: x[3], reverse=True)[:3]
        top3_anth = sorted(sens_list, key=lambda x: x[4], reverse=True)[:3]

        print(f'  FW Top 3: ' + ' | '.join([f'{n}({s:.3f})' for n,_,s,_,_ in top3_fw]))
        print(f'  Stress Top 3: ' + ' | '.join([f'{n}({s:.3f})' for n,_,_,s,_ in top3_st]))
        print(f'  Anth Top 3: ' + ' | '.join([f'{n}({s:.3f})' for n,_,_,_,s in top3_anth]))

    print()
    print('=' * 80)
    print('分析完成!')
    print('=' * 80)

    return all_results


if __name__ == '__main__':
    main()
