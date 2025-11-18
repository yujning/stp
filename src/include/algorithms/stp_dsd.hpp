#pragma once
#include <bits/stdc++.h>
using namespace std;

#include <chrono>
#include <fstream>
#include "excute.hpp"
#include"reorder.hpp"
void all_reorders(const std::string &binary);

//-----------------------------------------
// 生成单位矩阵表示
//-----------------------------------------
vector<stp_data> identity_vec(int dim)
{
    vector<stp_data> I(dim * dim + 1);
    I[0] = dim;
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            I[1 + i * dim + j] = (i == j) ? 1 : 0;
    return I;
}



//-----------------------------------------
// 打印矩阵(调试)
//-----------------------------------------
void print_vec(const vector<stp_data> &M)
{
    cout << "[";
    for (size_t i = 0; i < M.size(); ++i)
    {
        cout << M[i];
        if (i + 1 < M.size()) cout << " ";
    }
    cout << "]\n";
}





//-----------------------------------------
// 01 -> 12 表示（'0'->'2','1'->'1'）
//-----------------------------------------
static inline string to12(string b01){
    for(char& c: b01){ if(c=='0') c='2'; else c='1'; }
    return b01;
}

//-----------------------------------------
// 乘法 ui * w ： ui∈{ "12","21" }，w为任意偶长串
// [12]*w = w；[21]*w 每2位交换
//-----------------------------------------
static string mul_ui(const string& u, const string& w){
    // u ∈ {"12","21","11","22"}，w为偶数长度
    if(u=="12") return w;                 // 恒等
    if(u=="21"){                          // 每2位交换
        string r; r.reserve(w.size());
        for(size_t i=0;i<w.size(); i+=2){
            r.push_back(w[i+1]); r.push_back(w[i]);
        }
        return r;
    }
    if(u=="11" || u=="22"){               // 常量投影：忽略 w，直接返回全1/全2
        char c = (u=="11") ? '1' : '2';
        return string(w.size(), c);
    }
    return w;
}

//-----------------------------------------
// 输出每种情形说明
//-----------------------------------------
static string case_note(int cid){
    switch(cid){
        case 1: return "两个常量块（如全0与全1）";
        case 2: return "一种常量块 + 一种非常量块";
        case 3: return "只有一种非常量块";
        case 4: return "两种非常量块且互补";
        case 5: return "只有一种常量块（函数恒定）";
        default: return "未知";
    }
}

//==============================================================================
// 通用模板求解：按块分组→枚举 S → 构造 MF → 求/验 (MΦ, Mψ)
//==============================================================================

struct TemplateResult{
    string TypeName; // Type1/2/3/4
    string S0, S1;   // S = {S0, S1}
    string MF;       // [i1 i2 i3 i4] = [S0 S1]
    string Mphi;     // [j1 ... j_m]
    string Mpsi;     // ..
};

// blocks01: 原始01域分块；s：分块参数
static vector<TemplateResult> solve_with_template(const vector<string>& blocks01, int s)
{
    vector<string> W; W.reserve(blocks01.size());
    for(auto &b: blocks01) W.push_back(to12(b));
    const int m = (int)W.size();
    vector<TemplateResult> answers;

    // 按块值分组（最多两组）
    unordered_map<string,int> gid;
    vector<int> gidx(m, -1);
    vector<string> gval;
    for(int i=0;i<m;++i){
        if(!gid.count(W[i])){ gid[W[i]]=gval.size(); gval.push_back(W[i]); }
        gidx[i]=gid[W[i]];
    }
    if(gval.size()==0||gval.size()>2) return answers;

    bool has11=false, has22=false;
    for(auto &b: W){
        if(is_constant_block(b)){
            if(b[0]=='1') has11=true;
            if(b[0]=='2') has22=true;
        }
    }

    // 枚举所有可能的 S
    vector<pair<string,string>> Scands;
    Scands.push_back({"11","22"});
    Scands.push_back({"22","11"});
    if(has22){
        Scands.push_back({"22","12"});
        Scands.push_back({"22","21"});
        Scands.push_back({"21","22"});
        Scands.push_back({"12","22"});
    }
    if(has11){
        Scands.push_back({"11","12"});
        Scands.push_back({"11","21"});
        Scands.push_back({"12","11"});
        Scands.push_back({"21","11"});
    }
    Scands.push_back({"12","12"});
    Scands.push_back({"21","21"});
    Scands.push_back({"12","21"});
    Scands.push_back({"21","12"});

    auto is_invertible = [](const string& u){ return (u=="12"||u=="21"); };

    for(const auto& S : Scands)
    {
        string S0=S.first, S1=S.second;
        string MF=S0+S1;

        // 找两组代表
        int repA=-1,repB=-1;
        for(int i=0;i<m;++i){ if(gidx[i]==0){repA=i;break;} }
        for(int i=0;i<m;++i){ if(gidx[i]==1){repB=i;break;} }
        if(repA<0&&repB<0) continue;
        if(repA<0) repA=repB; if(repB<0) repB=repA;

        for(int map=0; map<2; ++map)
        {
            string u_g0=(map==0)?S0:S1;
            string u_g1=(map==0)?S1:S0;

            // === 修正版：只用可逆u反推Mψ，常量u只校验 ===
            vector<string> mpsi_candidates;
            if(is_invertible(u_g0)) mpsi_candidates.push_back(mul_ui(u_g0,W[repA]));
            if(is_invertible(u_g1)) mpsi_candidates.push_back(mul_ui(u_g1,W[repB]));

            string Mpsi;
            if(mpsi_candidates.empty()){
                int pick=-1;
                for(int i=0;i<m;++i) if(!is_constant_block(W[i])){pick=i;break;}
                if(pick<0) pick=repA;
                Mpsi=W[pick];
            }else{
                bool same=true;
                for(size_t k=1;k<mpsi_candidates.size();++k)
                    if(mpsi_candidates[k]!=mpsi_candidates[0]){same=false;break;}
                if(!same) continue;
                Mpsi=mpsi_candidates[0];
            }

            string gen0=mul_ui(u_g0,Mpsi);
            string gen1=mul_ui(u_g1,Mpsi);
            // ===================================================
            // 改进版：根据块分组与 j 分配动态构造 MΦ
            // ===================================================
            bool ok = true;
            string Mphi; Mphi.reserve(m);

            // 哪组取 j=1，哪组取 j=2
            // 默认 group0 对应 j=1, group1 对应 j=2，若反过来映射再交换
            for (int i = 0; i < m; ++i)
            {
                int gid_i = gidx[i];
                if ((gid_i == 0 && map == 0) || (gid_i == 1 && map == 1))
                    Mphi.push_back('1');
                else
                    Mphi.push_back('2');
            }

            // 逐块验证是否满足 Wi = [ip iq]*Mψ
            for (int i = 0; i < m; ++i)
            {
                string ipiq = (Mphi[i] == '1') ? S0 : S1;
                string expect = mul_ui(ipiq, Mpsi);
                if (W[i] != expect) { ok = false; break; }
            }
            if (!ok) continue;


            string typeName="Type 4";
            if((S0=="11"&&S1=="22")||(S0=="22"&&S1=="11")) typeName="Type 1";
            else if((S0=="12"&&S1=="12")||(S0=="21"&&S1=="21")) typeName="Type 3";
            else if((S0=="12"&&S1=="21")||(S0=="21"&&S1=="12")) typeName="Type 4";
            else typeName="Type 2";

            answers.push_back({typeName,S0,S1,MF,Mphi,Mpsi});
        }
    }
    return answers;
}
//-----------------------------------------
// 对给定 s 执行分块并仅在命中时打印 + 调用通用模板打印解
//-----------------------------------------
inline bool analyze_by_s(const std::string& binary, int s){
    const int n = std::log2(binary.size());
    const int r = n/2;
    if(!is_power_of_two(binary.size()) || (1<<n)!=(int)binary.size()) return false;
    if(s<1 || s>r) return false;

    const int block_len = 1<<s;
    const int num_blocks = binary.size() / block_len;
    std::vector<std::string> blocks01(num_blocks);
    for(int i=0;i<num_blocks;++i)
        blocks01[i]=binary.substr(i*block_len, block_len);

    int cid = theorem33_case_id(binary, s);
    if(cid==0) return false;

    cout << "\n===== 命中：s="<<s<<"，块长 2^s="<<block_len<<"，块数="<<num_blocks<<" =====\n";
    cout << "分块："; for(auto &b: blocks01) cout<<"["<<to12(b)<<"]";
    cout << "\n=> 情形("<<cid<<")：" << case_note(cid) << "\n";

    auto results = solve_with_template(blocks01, s);
    if(results.empty()){
        cout << "❌ 未找到可行的 (MF, MΦ, Mψ)\n";
    }else{
        int idx=1;
        for(auto &r:results){
            cout << "\n【"<<r.TypeName<<"】 S = {δ2["<<r.S0<<"], δ2["<<r.S1<<"]}\n";
            cout << idx++ << ". MF = ["<<r.MF<<"]\n";
            cout << "   MΦ = ["<<r.Mphi<<"]\n\n";
            cout << "   Mψ = ["<<r.Mpsi<<"]\n";
        }
    }
    return true;
}
//-----------------------------------------
// 批量检测（只打印符合的）
//-----------------------------------------
void analyze_all_s(const std::string& binary){
    const int n = std::log2(binary.size());
    if(!is_power_of_two(binary.size()) || (1<<n)!=(int)binary.size()){
        std::cout << "输入长度必须是 2 的整数次幂\n"; return;
    }
    const int r = n/2;
    for(int s=1; s<=r; ++s)
        analyze_by_s(binary, s);
}

//-----------------------------------------
// 主过程：生成所有重排并分析（仅打印符合的）
//-----------------------------------------
