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


// =============================
// 通用模板分解
// =============================
struct TemplateResult {
    string MF;
    string Mphi;
    string Mpsi;
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

    // ----------------- 求 MΦ -----------------
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
// 创建节点
// =============================
static int new_node(const string& func, const vector<int>& child)
{
    int id = NODE_ID++;
    NODE_LIST.push_back({id, func, child});
    return id;
}


// =============================
// 构造小 LUT 的“叶子结构”
// 长度 <=4 时不直接当叶，而是再挂一层输入节点：
//  len==2:  F -> F(var)
//  len==4:  F -> F(varL, varR)
// =============================
static int build_small_tree(const string& f12)
{
    int len = f12.size();

    if(len == 2){
        // 建一个输入节点（真值表随便，标个 "in" 就行）
        int in = new_node("in", {});         // 最底层叶子
        int op = new_node(f12, {in});        // 一元结点，如 12(in) 或 21(in)
        return op;
    }
    else if(len == 4){
        int inL = new_node("in", {});        // 左输入
        int inR = new_node("in", {});        // 右输入
        int op  = new_node(f12, {inL, inR}); // 二元结点，如 1222(inL, inR)
        return op;
    }
    else{
        // 理论上不会来这里，保险起见
        return new_node(f12, {});
    }
}


// =============================
// 递归分解（长度≤4 用 build_small_tree）
// =============================
static int decompose(const string& f12, int depth=0)
{
    int len = f12.size();

    // ----------- 小 LUT：长度 ≤ 4 -----------
    if(len <= 4){
        // 不再往下做定理 3.3 分解，直接挂输入叶节点
        return build_small_tree(f12);
    }

    string f01 = to01(f12);
    if(!is_power_of_two(f01.size())){
        return build_small_tree(f12);
    }

    int n = log2(len);
    int r = n/2;

    // 尝试所有 s
    for(int s=1;s<=r;s++){
        int cid = theorem33_case_id(f01, s);
        if(cid==0) continue;

        int bl = 1<<s;
        int nb = len/bl;

        vector<string> blocks01(nb);
        for(int i=0;i<nb;i++)
            blocks01[i] = f01.substr(i*bl, bl);

        // 判断常量类型
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
                S_list={{"11","22"},{"22","11"}};
                break;
            case 2:
                if(has11)
                    S_list={{"11","12"},{"11","21"},{"12","11"},{"21","11"}};
                else
                    S_list={{"22","12"},{"22","21"},{"12","22"},{"21","22"}};
                break;
            case 3:
                S_list={{"12","12"},{"21","21"}};
                break;
            case 4:
                S_list={{"12","21"},{"21","12"}};
                break;
            case 5:
                return build_small_tree(f12);
        }

        // 只取第一个 S（按你的 case 规则）
        auto ret = run_case_once(blocks01, s, S_list[0].first, S_list[0].second);

        // ------ 打印此层分解结果 ------
        indent(depth);
        cout << "MF  = [" << ret.MF  << "]\n";
        indent(depth);
        cout << "MΦ  = [" << ret.Mphi << "]\n";
        indent(depth);
        cout << "Mψ  = [" << ret.Mpsi << "]\n\n";

        // ------ 递归（Mphi/Mpsi 里如果长度≤4，会走 build_small_tree）------
        int left_id  = decompose(ret.Mphi, depth+1);
        int right_id = decompose(ret.Mpsi, depth+1);

        // ------ 构造当前节点 ------
        return new_node(ret.MF, {left_id, right_id});
    }

    // 无法分解 → 当小 LUT 处理
    return build_small_tree(f12);
}


// =============================
// 顶层入口：analyze_by_s
// =============================
inline bool analyze_by_s(const string& binary, int s)
{
    int len = binary.size();
    if(!is_power_of_two(len)) return false;

    int bl = 1<<s;
    int nb = len/bl;

    vector<string> blocks01(nb);
    for(int i=0;i<nb;i++)
        blocks01[i]=binary.substr(i*bl, bl);

    int cid = theorem33_case_id(binary, s);
    if(cid==0) return false;

    cout << "\n===== 命中 s="<<s<<" =====\n";
    cout << "分块：";
    for(auto& b:blocks01) cout<<"["<<to12(b)<<"]";
    cout<<"\n\n";

    // 重置节点表
    NODE_LIST.clear();
    NODE_ID = 1;

    string f12 = to12(binary);

    // ------ 构造 root ------
    int root = decompose(f12, 0);

    // ------ 打印所有节点 ------
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
