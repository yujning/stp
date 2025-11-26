#pragma once
#include <bits/stdc++.h>
using namespace std;

#include "excute.hpp"
#include "reorder.hpp"

// =============================
// DSD 节点结构
// =============================
struct DSDNode {
    int id;                 // 节点编号
    string func;            // 12 域真值表
    vector<int> child;      // 子节点编号
};

static vector<DSDNode> NODE_LIST;
static int NODE_ID = 1;
static int STEP_ID = 1;     // 分解步骤编号：1,2,3,...

// =============================
// 工具函数
// =============================
static inline string to12(const string& b01)
{
    string r=b01;
    for(char& c:r)
        c = (c=='0' ? '2' : '1');
    return r;
}

static inline string to01(const string& b12)
{
    string r=b12;
    for(char& c:r)
        c = (c=='1' ? '1' : '0');
    return r;
}

// ui ∈ {11,12,21,22}，w 偶长
static inline string mul_ui(const string& ui, const string& w)
{
    if(ui=="12") return w;              // 恒等
    if(ui=="21"){                       // swap
        string r; r.reserve(w.size());
        for(size_t i=0;i<w.size();i+=2){
            r.push_back(w[i+1]);
            r.push_back(w[i]);
        }
        return r;
    }
    if(ui=="11" || ui=="22"){           // 常量投影
        char c = (ui=="11" ? '1' : '2');
        return string(w.size(), c);
    }
    return w;
}

static void indent(int d)
{
    while(d--) cout<<"   ";
}

static string case_note(int cid)
{
    switch(cid){
        case 1: return "两个常量块（如全0与全1）";
        case 2: return "一种常量块 + 一种非常量块";
        case 3: return "只有一种非常量块";
        case 4: return "两种非常量块且互补";
        case 5: return "只有一种常量块（函数恒定）";
        default: return "未知";
    }
}

// =============================
// 通用模板分解：给定 blocks01 + S0,S1
// =============================
struct TemplateResult {
    string MF;
    string Mphi;   // 标签串：每位 '1' 或 '2'
    string Mpsi;   // 真实 12 域真值表
};

static TemplateResult run_case_once(
    const vector<string>& blocks01,
    int s,
    const string& S0,
    const string& S1)
{
    vector<string> W;
    for(auto& b:blocks01)
        W.push_back(to12(b));

    int m = W.size();
    string MF = S0 + S1;

    // ----------------- 求 Mψ -----------------
    string Mpsi;
    bool found=false;

    auto try_u=[&](const string& u)->bool{
        if(u=="12" || u=="21"){
            vector<string> cand;

            for(int i=0;i<m;i++){
                // 跳过常量块
                if(W[i] == string(W[i].size(), W[i][0]))
                    continue;

                string c = mul_ui(u, W[i]);

                // 检查是否 c 反推回去能得到 W[i]
                if(mul_ui(u, c) == W[i])
                    cand.push_back(c);
            }

            if(cand.empty()) return false;

            // 多个非常量块反推结果必须一致
            for(size_t k=1;k<cand.size();k++)
                if(cand[k] != cand[0])
                    return false;

            Mpsi = cand[0];
            return true;
        }
        return false;
    };

    if(!found && try_u(S0)) found=true;
    if(!found && try_u(S1)) found=true;

    if(!found){
        int pick=-1;
        for(int i=0;i<m;i++){
            if(!is_constant_block(blocks01[i])){
                pick=i; break;
            }
        }
        if(pick<0) pick=0;
        Mpsi = W[pick];
    }

    // ----------------- 求 MΦ 标签 -----------------
    string exp0 = mul_ui(S0, Mpsi);
    string exp1 = mul_ui(S1, Mpsi);
    string Mphi; Mphi.reserve(m);

    for(int i=0;i<m;i++){
        if(W[i] == exp0) Mphi.push_back('1');
        else             Mphi.push_back('2');
    }

    return {MF, Mphi, Mpsi};
}

// =============================
// 创建节点 + 小 LUT 叶子结构
// =============================
static int new_node(const string& func, const vector<int>& child)
{
    int id = NODE_ID++;
    NODE_LIST.push_back({id, func, child});
    return id;
}

// 长度 ≤4：构造 2-LUT 叶子结构
static int build_small_tree_from01(const string& f01)
{
    string f12 = to12(f01);
    int len = f12.size();

    if(len == 2){
        int in = new_node("in", {});
        return new_node(f12, {in});
    }
    else if(len == 4){
        int inL = new_node("in", {});
        int inR = new_node("in", {});
        return new_node(f12, {inL, inR});
    }
    // 理论上不会进来
    return new_node(f12, {});
}

// =============================
// 单步：对一个 01 串做 “重排 + case + 模板”
// 输入：binary01（当前要分解的函数，01 域）
// 输出：MF12（4位12串），phi01, psi01（作为下一轮的 01 域真值表）
// 返回：是否成功分解（false = 无法按定理3.3分解）
// =============================
static bool factor_once_with_reorder_01(
    const string& binary01,
    int depth,
    string& MF12,
    string& phi01,
    string& psi01)
{
    int len = binary01.size();
    if(!is_power_of_two(len) || len <= 4) return false;

    int n = log2(len);
    int r = n/2;

    vector<stp_data> Mf = binary_to_vec(binary01);

    for(int s=1; s<=r; ++s)
    {
        vector<bool> v(n);
        fill(v.begin(), v.begin()+s, true);

        do{
            // Λ
            vector<int> Lambda;
            for(int i=0;i<n;i++)
                if(v[i]) Lambda.push_back(i+1);

            // swap_chain
            vector<vector<stp_data>> swap_chain;
            for(int k=s;k>=1;k--){
                int j_k = Lambda[k-1];
                int exp = j_k + (s-1) - k;
                swap_chain.push_back(generate_swap_vec(2, pow(2,exp)));
            }
            swap_chain.push_back(generate_swap_vec(pow(2,n-s), pow(2,s)));

            vector<stp_data> Mperm =
                Vec_chain_multiply(swap_chain,false);

            vector<stp_data> result =
                Vec_semi_tensor_product(Mf,Mperm);

            string reordered;
            reordered.reserve(len);
            for(size_t i=1;i<result.size();++i)
                reordered.push_back(result[i]?'1':'0');

            int cid = theorem33_case_id(reordered, s);
            if(cid!=0){
                // ---------- 打印重排命中 ----------
                // cout << "\n===== 重排命中：s="<<s<<" 情形("<<cid<<") =====\n";
                // cout << "Λ = { ";
                // for(int j : Lambda) cout<<j<<" ";
                // cout << "}  => reordered: " << reordered << "\n";

                // ---------- 分块 ----------
                int bl = 1<<s;
                int nb = len/bl;
                vector<string> blocks01(nb);
                for(int i=0;i<nb;i++)
                    blocks01[i] = reordered.substr(i*bl, bl);

               // cout << "\n===== 命中：s="<<s<<"，块长="<<bl<<"，块数="<<nb<<" =====\n";
                //cout << "分块：";
                //for(auto& b:blocks01) cout<<"["<<to12(b)<<"]";
                //cout << "\n=> 情形("<<cid<<")："<<case_note(cid)<<"\n\n";

                // ---------- 按 case 选 S0,S1 ----------
                bool has11=false, has22=false;
                for(auto& b:blocks01){
                    if(is_constant_block(b)){
                        if(b[0]=='1') has11=true;
                        if(b[0]=='0') has22=true;
                    }
                }

                vector<pair<string,string>> S_list;
                switch(cid){
                    case 1:
                        S_list = {{"11","22"},{"22","11"}};
                        break;
                    case 2:
                        if(has11)
                            S_list = {{"11","12"},{"11","21"},{"12","11"},{"21","11"}};
                        else
                            S_list = {{"22","12"},{"22","21"},{"12","22"},{"21","22"}};
                        break;
                    case 3:
                        S_list = {{"12","12"},{"21","21"}};
                        break;
                    case 4:
                        S_list = {{"12","21"},{"21","12"}};
                        break;
                    case 5:
                        return false;   // 恒定函数，不再分解
                }

                string S0 = S_list[0].first;
                string S1 = S_list[0].second;

                auto R = run_case_once(blocks01, s, S0, S1);

                // ---------- 打印这一步分解 ----------
                cout << STEP_ID++ << ". MF  = [" << R.MF  << "]\n";
                cout << "   MΦ  = [" << R.Mphi << "]\n";
                cout << "   Mψ  = [" << R.Mpsi << "]\n\n";

                MF12 = R.MF;
                phi01 = to01(R.Mphi);   // 标签串 1/2 -> 01
                psi01 = to01(R.Mpsi);   // 真值表 12 -> 01
                return true;
            }

        }while(prev_permutation(v.begin(), v.end()));
    }

    return false;
}

// =============================
// 主递归：对 01 域真值表做 DSD
// =============================
static int dsd_factor(const string& f01, int depth=0)
{
    int len = f01.size();

    // 小函数：直接构造 2-LUT 结构
    if(len <= 4 || !is_power_of_two(len)){
        return build_small_tree_from01(f01);
    }

    string MF12, phi01, psi01;
    if(!factor_once_with_reorder_01(f01, depth, MF12, phi01, psi01)){
        // 不能按定理 3.3 分解，当成 2-LUT
        return build_small_tree_from01(f01);
    }

    int left_id  = dsd_factor(phi01, depth+1);
    int right_id = dsd_factor(psi01, depth+1);

    return new_node(MF12, {left_id, right_id});
}

// =============================
// 顶层入口：递归 DSD
// =============================
inline bool run_dsd_recursive(const string& binary01)
{
    if(!is_power_of_two(binary01.size())){
        cout<<"输入长度必须是 2 的整数次幂\n";
        return false;
    }

    NODE_LIST.clear();
    NODE_ID = 1;
    STEP_ID = 1;

    cout << "输入 = " << binary01
         << " (len = " << binary01.size() << ")\n";

    int root = dsd_factor(binary01, 0);

    cout << "===== 最终 DSD 节点列表 =====\n";
    for(auto& nd : NODE_LIST){
        cout << nd.id << " = " << nd.func;
        if(nd.child.size()==1){
            cout<<"("<<nd.child[0]<<")";
        }else if(nd.child.size()==2){
            cout<<"("<<nd.child[0]<<","<<nd.child[1]<<")";
        }
        cout<<"\n";
    }
    cout << "Root = " << root << "\n";
    return true;
}

// =============================
// 兼容接口：给 reorder.hpp 里 all_reorders 调用
// （如果还在用的话）
// 这里简单地把 binary 当新的 f01 丢给 run_dsd_recursive
// s 已经没用了，仅保持签名不报错
// =============================
inline bool analyze_by_s(const string& binary, int s)
{
    (void)s;
    return run_dsd_recursive(binary);
}
