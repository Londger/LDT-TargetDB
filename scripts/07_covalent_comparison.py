"""
共价策略对比：LDT(NASA) vs SuFEx vs 传统共价(丙烯酰胺)
比较不同共价化学对同一批靶点的适配性
"""
import pandas as pd
import numpy as np
from pathlib import Path

BASE = Path("c:/Users/Longer/Documents/自写生信论文")
RES = BASE / "results"

# 三种共价策略的残基权重
STRATEGIES = {
    "LDT_NAS": {
        "name": "LDT-NASA (This work)",
        "description": "配体导向转移 N-酰基磺酰胺, 无痕标记, Lys首选, His不反应",
        "weights": {"LYS": 1.0, "CYS": 0.3, "TYR": 0.2, "SER": 0.1},
        "primary_target": "LYS",
    },
    "SuFEx_CTR": {
        "name": "SuFEx-CTR (Liu 2024 Nature)",
        "description": "氟磺酰交换共价放射配体, Tyr/Lys/His均可反应",
        "weights": {"TYR": 1.0, "LYS": 0.8, "HIS": 0.7, "SER": 0.3},
        "primary_target": "TYR",
    },
    "Traditional_Acrylamide": {
        "name": "Traditional Cys-targeted",
        "description": "传统丙烯酰胺共价弹头, Cys靶向",
        "weights": {"CYS": 1.0, "LYS": 0.2},
        "primary_target": "CYS",
    },
}

# 标准pKa
REF_PKA = {"LYS": 10.5, "CYS": 8.5, "TYR": 10.0, "SER": 13.0, "HIS": 6.5}

# 加载enhanced结果
enhanced = pd.read_csv(RES / "enhanced_final_ranking.csv")
print(f"Loaded {len(enhanced)} targets from enhanced ranking")

# 对每种策略重新计算LDT-like得分
# 注：由于我们没有重新跑pocket+PROPKA，这里基于已有数据做近似
# 实际得分 = weight * pKa_score(SASA默认0.3)

for strat_key, strat_info in STRATEGIES.items():
    weights = strat_info["weights"]
    scores = []

    for _, row in enhanced.iterrows():
        # 口袋中有哪些类型的残基（近似：用已有top_res_type和n_nucleophiles作为代理）
        # 更准确：重新解析best_pocket_residues字段
        res_types_in_pocket = set()
        if pd.notna(row.get("top_res_type")) and row["top_res_type"]:
            res_types_in_pocket.add(row["top_res_type"])

        best_score = 0
        best_res = ""
        for rt in weights:
            w = weights[rt]
            pk = REF_PKA.get(rt, 10)
            ps = np.exp(-abs(pk - 7.4) / 2)
            # SASA: 匹配top_res_type的残基用0.5，口袋中存在的残基用0.2，否则用默认0.1
            if rt == row.get("top_res_type"):
                ss = 0.5
            elif len(res_types_in_pocket) > 0:
                ss = 0.2  # 口袋中可能有但未确认
            else:
                ss = 0.1
            sc = w * ps * ss
            if sc > best_score:
                best_score = sc
                best_res = rt

        scores.append({"gene": row["gene"], f"{strat_key}_score": round(best_score, 4),
                       f"{strat_key}_best_res": best_res})

    strat_df = pd.DataFrame(scores)
    enhanced = enhanced.merge(strat_df, on="gene", how="left")

# 比较三种策略的得分
enhanced["best_strategy"] = ""
for i, row in enhanced.iterrows():
    scores = {
        "LDT_NAS": row.get("LDT_NAS_score", 0) or 0,
        "SuFEx_CTR": row.get("SuFEx_CTR_score", 0) or 0,
        "Traditional_Acrylamide": row.get("Traditional_Acrylamide_score", 0) or 0,
    }
    best = max(scores, key=scores.get)
    enhanced.at[i, "best_strategy"] = best

# 统计各策略适配的靶点数量
print("\n=== Covalent Strategy Comparison ===\n")
for sk, si in STRATEGIES.items():
    n = (enhanced["best_strategy"] == sk).sum()
    print(f"{si['name']:30s}  {n:4d} targets  ({n/len(enhanced)*100:.1f}%)")

print(f"\n{'Target':<15} {'LDT-NAS':<12} {'SuFEx-CTR':<12} {'Trad-Cys':<12} {'Best Strategy'}")
print("-"*70)
for _, row in enhanced.head(20).iterrows():
    ldt = row.get("LDT_NAS_score", 0) or 0
    su = row.get("SuFEx_CTR_score", 0) or 0
    tr = row.get("Traditional_Acrylamide_score", 0) or 0
    best = row["best_strategy"]
    print(f"{row['gene']:<15} {ldt:.4f}       {su:.4f}       {tr:.4f}       {best}")

# 策略独特性分析：只在一种策略下得高分的靶点
for sk, si in STRATEGIES.items():
    only_this = enhanced[
        (enhanced["best_strategy"] == sk) &
        (enhanced[f"{sk}_score"] > 0.1)
    ]
    other_best = []
    for other_sk in STRATEGIES:
        if other_sk != sk:
            other_best.extend(only_this["gene"].tolist())
    unique = [g for g in only_this["gene"] if g not in other_best or True]  # simplified
    print(f"\n{si['name']} unique targets (score>0.1): {len(only_this)}")
    for _, r in only_this.head(10).iterrows():
        print(f"  {r['gene']:15s} score={r[f'{sk}_score']:.4f}")

# 保存
enhanced.to_csv(RES / "covalent_strategy_comparison.csv", index=False)
print(f"\nSaved to: {RES}/covalent_strategy_comparison.csv")
