#pragma once
#include <bits/stdc++.h>
using namespace std;

#include <chrono>
#include <fstream>
#include "excute.hpp"
#include "reorder.hpp"   // is_constant_block / theorem33_case_id / all_reorders
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
// 12 -> 01 表示（'1'->'1','2'->'0'）
// 用于把 MΦ / Mψ 当成新的 01 域函数再分解
//-----------------------------------------
static inline string to01_from12(string w12){
    string r = w12;
    for(char& c : r) c = (c=='1' ? '1' : '0');
    return r;
}

//-----------------------------------------
// 乘法 ui * w ： ui∈{ "12","21","11","22" }，w为任意偶长串
// [12]*w = w；[21]*w 每2位交换；[11]/[22] 常量投影
//-----------------------------------------
static string mul_ui(const string& u, const string& w){
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

//-----------------------------------------
// 模板结果结构体
//-----------------------------------------
struct TemplateResult{
    string S0, S1;   // 输入的 S = {S0, S1}
    string MF;       // S0S1 拼接
    string Mphi;     // 1/2 字符串（12 域）
    string Mpsi;     // 1/2 字符串（12 域）
};

//-----------------------------------------
// 固定 S 情况下的通用模板求解
// blocks01 是 "01 域"，内部会转成 12 域
//-----------------------------------------
static TemplateResult run_template_with_fixed_S(
    const vector<string>& blocks01,
    int /*s*/,
    const string& S0,
    const string& S1
){
    vector<string> W;
    for (auto& b : blocks01) W.push_back(to12(b));
    int m = W.size();

    string MF = S0 + S1;

    // 1. 找 Mψ
    string Mpsi;
    {
        bool found = false;
        auto try_u = [&](const string& u){
            if (u=="12" || u=="21"){
                for(int i=0;i<m;i++){
                    string candidate = mul_ui(u, W[i]);
                    if (mul_ui(u, candidate) == W[i]) {
                        Mpsi = candidate;
                        return true;
                    }
                }
            }
            return false;
        };

        if (!found && try_u(S0)) found=true;
        if (!found && try_u(S1)) found=true;

        // 若还是没找到，取第一个非常量块
        if (!found){
            int pick=-1;
            for(int i=0;i<m;i++){
                if (!is_constant_block(blocks01[i])){
                    pick=i; break;
                }
            }
            if (pick<0) pick=0;
            Mpsi = W[pick];
        }
    }

    // 2. 求 MΦ
    string Mphi; Mphi.reserve(m);
    string expect0 = mul_ui(S0, Mpsi);
    string expect1 = mul_ui(S1, Mpsi);
    for (int i=0;i<m;i++){
        if (W[i] == expect0) Mphi.push_back('1');
        else                 Mphi.push_back('2');
    }

    return { S0, S1, MF, Mphi, Mpsi };
}

//-----------------------------------------
// 打印缩进（递归用）
//-----------------------------------------
static inline void print_indent(int depth){
    for(int i=0;i<depth;i++) cout << "   ";
}

//-----------------------------------------
// 递归分解：对 12 域函数 f12 继续分解到 2-LUT
//-----------------------------------------
static void recursive_decompose_12(const string& f12, int depth)
{
    string bin = to01_from12(f12);
    int len = bin.size();

    // 终止条件：长度 <= 4 -> 2-LUT 叶子
    if(len <= 4){
        print_indent(depth);
        cout << "Leaf (2-LUT): f = [" << f12 << "] (len=" << len << ")\n";
        return;
    }

    if(!is_power_of_two(len)){
        print_indent(depth);
        cout << "无法继续分解（长度非 2 的幂）: f = [" << f12 << "] len=" << len << "\n";
        return;
    }

    int n = std::log2(len);
    int r = n / 2;
    bool found = false;

    for(int s=1; s<=r && !found; ++s)
    {
        int cid = theorem33_case_id(bin, s);
        if (cid == 0) continue;

        int block_len = 1 << s;
        int num_blocks = len / block_len;
        vector<string> blocks01(num_blocks);
        for(int i=0;i<num_blocks;i++)
            blocks01[i] = bin.substr(i*block_len, block_len);

        // 统计常量类型
        bool has11=false, has22=false;
        for(auto& b : blocks01){
            if(is_constant_block(b)){
                if (b[0]=='1') has11=true; // 全1块
                if (b[0]=='0') has22=true; // 全0块（映射到 22）
            }
        }

        vector<pair<string,string>> S_list;

        // --- 与顶层 analyze_by_s 一致的 Case→S 映射 ---
        if (cid == 1)
        {
            S_list = { {"11","22"}, {"22","11"} };
        }
        else if (cid == 2)
        {
            if (has11){ // 常量 = 11
                S_list = {
                    {"11","12"}, {"11","21"},
                    {"12","11"}, {"21","11"}
                };
            }else{ // 常量 = 22
                S_list = {
                    {"22","12"}, {"22","21"},
                    {"12","22"}, {"21","22"}
                };
            }
        }
        else if (cid == 3)
        {
            S_list = { {"12","12"}, {"21","21"} };
        }
        else if (cid == 4)
        {
            S_list = { {"12","21"}, {"21","12"} };
        }
        else if (cid == 5)
        {
            print_indent(depth);
            cout << "函数恒定（Case 5），停止分解: f = ["<< f12 <<"]\n";
            return;
        }

        // 这里我们取第一个可行的 S（你以后也可以改成枚举所有 S）
        for(auto& S : S_list)
        {
            TemplateResult R = run_template_with_fixed_S(blocks01, s, S.first, S.second);

            print_indent(depth);
            cout << "递归分解：s="<<s<<" case="<<cid
                 << "  MF=["<<R.MF<<"]  对应 f=["<<f12<<"]\n";
            print_indent(depth);
            cout << "   MΦ = [" << R.Mphi << "]\n";
            print_indent(depth);
            cout << "   Mψ = [" << R.Mpsi << "]\n";

            // 递归处理 MΦ
            if(to01_from12(R.Mphi).size() > 4){
                recursive_decompose_12(R.Mphi, depth+1);
            }else{
                print_indent(depth+1);
                cout << "Φ 叶子 (2-LUT): ["<< R.Mphi <<"]\n";
            }

            // 递归处理 Mψ
            if(to01_from12(R.Mpsi).size() > 4){
                recursive_decompose_12(R.Mpsi, depth+1);
            }else{
                print_indent(depth+1);
                cout << "ψ 叶子 (2-LUT): ["<< R.Mpsi <<"]\n";
            }

            found = true;
            break;
        }
    }

    if(!found){
        print_indent(depth);
        cout << "未找到可用的 Case 继续分解: f = ["<< f12 <<"]\n";
    }
}

//-----------------------------------------
// 顶层：对给定 s 执行分块 + Case 判定 + 固定 S + 通用模板
// 然后对 MΦ/Mψ 做递归分解
//-----------------------------------------
inline bool analyze_by_s(const std::string& binary, int s){
    const int n = std::log2(binary.size());
    if(!is_power_of_two(binary.size()) || (1<<n)!=(int)binary.size()) return false;

    int r = n/2;
    if (s < 1 || s > r) return false;

    int block_len = 1 << s;
    int num_blocks = binary.size() / block_len;

    vector<string> blocks01(num_blocks);
    for(int i=0;i<num_blocks;i++)
        blocks01[i] = binary.substr(i*block_len, block_len);

    int cid = theorem33_case_id(binary, s);
    if (cid == 0) return false;

    cout << "\n===== 命中：s="<<s<<"，块长="<<block_len<<"，块数="<<num_blocks<<" =====\n";
    cout << "分块：";
    for(auto& b: blocks01) cout<<"["<<to12(b)<<"]";
    cout << "\n=> 情形("<<cid<<")：" << case_note(cid) << "\n\n";

    vector<pair<string,string>> S_list;

    // 统计常量类型（注意 01 域：0 -> 22, 1 -> 11）
    bool has11=false, has22=false;
    for(auto& b : blocks01){
        if(is_constant_block(b)){
            if (b[0]=='1') has11=true; // 全1块 -> 11
            if (b[0]=='0') has22=true; // 全0块 -> 22
        }
    }

    // =============================
    // Case 1：两个常量块
    // =============================
    if (cid == 1)
    {
        S_list = { {"11","22"}, {"22","11"} };
    }

    // =============================
    // Case 2：一种常量 + 一种非常量
    // =============================
    else if (cid == 2)
    {
        if (has11){ // 常量 = 11
            S_list = {
                {"11","12"}, {"11","21"},
                {"12","11"}, {"21","11"}
            };
        }else{ // 常量 = 22
            S_list = {
                {"22","12"}, {"22","21"},
                {"12","22"}, {"21","22"}
            };
        }
    }

    // =============================
    // Case 3：只有一种非常量
    // =============================
    else if (cid == 3)
    {
        S_list = { {"12","12"}, {"21","21"} };
    }

    // =============================
    // Case 4：两个互补非常量
    // =============================
    else if (cid == 4)
    {
        S_list = { {"12","21"}, {"21","12"} };
    }

    // =============================
    // Case 5：恒定函数
    // =============================
    else if (cid == 5)
    {
        cout << "函数恒定，无意义。\n";
        return true;
    }

    // =============================
    // 进行计算：使用每个 S
    // =============================
    int idx=1;
    for(auto& S : S_list)
    {
        auto R = run_template_with_fixed_S(blocks01, s, S.first, S.second);

        cout << idx++ << ". 顶层 MF = [" << R.MF << "]\n";
        cout << "   MΦ = [" << R.Mphi << "]\n";
        cout << "   Mψ = [" << R.Mpsi << "]\n\n";

        // === 递归分解 MΦ / Mψ 到 2-LUT ===
        if(to01_from12(R.Mphi).size() > 4 || to01_from12(R.Mpsi).size() > 4){
            cout << "   ⇩ 对子函数继续分解：\n";
            if(to01_from12(R.Mphi).size() > 4)
                recursive_decompose_12(R.Mphi, 2); // depth=2，缩进更深一点
            else{
                print_indent(2);
                cout << "Φ 已是 2-LUT: ["<< R.Mphi <<"]\n";
            }

            if(to01_from12(R.Mpsi).size() > 4)
                recursive_decompose_12(R.Mpsi, 2);
            else{
                print_indent(2);
                cout << "ψ 已是 2-LUT: ["<< R.Mpsi <<"]\n";
            }
        }
    }

    return true;
}

//-----------------------------------------
// 批量检测（只打印符合的）
//-----------------------------------------
inline void analyze_all_s(const std::string& binary){
    const int n = std::log2(binary.size());
    if(!is_power_of_two(binary.size()) || (1<<n)!=(int)binary.size()){
        std::cout << "输入长度必须是 2 的整数次幂\n"; return;
    }
    const int r = n/2;
    for(int s=1; s<=r; ++s)
        analyze_by_s(binary, s);
}
