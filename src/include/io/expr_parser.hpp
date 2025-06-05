#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>
#include <io/excute.hpp>


namespace stp
{
    enum NODE_TYPE
    {
        NodeType_None,
        NodeType_Gate,
        NodeType_Variable
    };

    enum GATE_TYPE
    {
        GateType_NONE,
        GateType_W,
        GateType_I,
        GateType_Mr,
        GateType_Mc,
        GateType_Md,
        GateType_Mn
    };

    class expr_node
    {
    public:
        //init variable node
        expr_node(const NODE_TYPE &nt , const std::string &name, const stp_data &id, const stp_data &k = 2)
            : nodetype(nt), node_name(name), var_id(id), var_k(k) {}

        //init gate node
        expr_node(const NODE_TYPE &nt, const std::string &name, const GATE_TYPE &gt, const stp_data &dim1, const stp_data &dim2)
            : nodetype(nt), gatetype(gt), node_name(name), dim_1(dim1), dim_2(dim2) {}

        // get node type
        const NODE_TYPE &Get_NodeType() const { return nodetype; }
        // get gate type
        const GATE_TYPE &Get_GateType() const { return gatetype; }
        // get variable id
        id Get_Var_id() const { return var_id; }
        // get variable k
        id Get_Var_k() const { return var_k; }
        //get dimension 1
        stp_data Get_dim1() const { return dim_1; }
        // get parameter 1
        stp_data Get_dim2() const { return dim_2; }
        // get node name
        std::string Get_Node_Name() const { return node_name; }

        void Set_dim2(stp_data d2) { dim_2 = d2; }

    private:
        NODE_TYPE nodetype;
        GATE_TYPE gatetype;
        std::string node_name;
        stp_data dim_1;
        stp_data dim_2;
        id var_id;
        stp_data var_k;
    };

    class expr_parser
    {
    public:
        //initialize
        expr_parser(const std::vector<std::string> &input_names, const std::vector<stp_data> &input_k, const std::string &expression, const bool& verbose = false)
                    : input_names(input_names),expression(expression),input_k(input_k),verbose(verbose)
        {
            //record the order of input variables
            input_to_id();
            from_exp_to_nmx();
            normalize_expr_chain();
            from_expr_to_matrix();
            exchange_vars_mixed();
        }

        bool Is_Variable(const expr_node n) { return (n.Get_NodeType() == NodeType_Variable); }
        bool Is_Gate(const expr_node n) { return (n.Get_NodeType() == NodeType_Gate); }
        GATE_TYPE Get_GateType(const expr_node n) { return n.Get_GateType(); }
        std::vector<expr_node> Get_Expr_Chain() const { return expr_chain; }

        void input_to_id(void)
        {
            for (const auto &name : this->input_names)
            {
                auto it = nameToId.find(name);
                //existed
                if (it != nameToId.end())
                {
                    //error
                    std::cout << "error: variable " << name << " already exists" << std::endl;
                }
                else
                {
                    id index = nameToId.size();
                    nameToId[name] = index;
                }
            }
        }

        // example:     (a & b) | (a & ~c) | (~b & ~c)   a b c
        void from_exp_to_nmx(void)
        {
            std::vector<std::string> equation;
            equation.push_back("(");
            std::string temp{""};

            for (int i = 0; i < expression.size(); i++)
            {
                if (expression[i] == ' ')
                {
                    // continue;
                }
                else if (expression[i] == '(' || expression[i] == ')' ||
                         expression[i] == '|' || expression[i] == '&' ||
                         expression[i] == '~')
                {
                    if (!temp.empty())
                    {
                        equation.push_back(temp);
                    }

                    std::string p1 = "";
                    p1 = p1 + expression[i];
                    equation.push_back(p1);
                    temp = "";
                }
                else
                {
                    temp = temp + expression[i];
                }
            }
            equation.push_back(")");

            // equation = ["(",   "(", "a", "&", "b", ")",   "(", "a", "&", "~", "c", ")",
            // "(", "~", "a", "&", "~", "c", ")"   ")"]

            /***********/
            // deal with "~"
            for (int i = 0; i < equation.size(); i++)
            {
                if (equation[i] == "&" || equation[i] == "|" ||
                    equation[i] == "(" || equation[i] == ")")
                {
                    // continue;
                }
                else if (equation[i] == "~")
                {
                    equation[i] = "~(" + equation[i + 1] + ")";
                    equation.erase(equation.begin() + i + 1);
                }
                else
                {
                    equation[i] = "(" + equation[i] + ")";
                }
            }

            // equation = ["(",   "(", "(a)", "&", "(b)", ")",   "(", "(a)", "&", "~(c)",
            // ")",  "(", "~(a)", "&", "~(c)", ")"   ")"]
            std::vector<int> left_bracket;
            for (int i = 0; i < equation.size(); i++)
            {
                if (equation[i] == "(")
                {
                    left_bracket.push_back(i);
                }
                if (equation[i] == ")")
                {
                    std::string equ = "";
                    for (int j = left_bracket[left_bracket.size() - 1] + 1; j < i;
                         j++) // 遍历左右括号之间的表达式
                    {
                        if (equation[j] == "|" || equation[j] == "&")
                        {
                            equ = "(" + equation[j] + ")" + equ;
                        }
                        else
                        {
                            equ = equ + "(" + equation[j] + ")";
                        }
                    }

                    //"(" -> equ
                    equation[left_bracket[left_bracket.size() - 1]] = equ;
                    equation.erase(
                        equation.begin() + left_bracket[left_bracket.size() - 1] + 1,
                        equation.begin() + i + 1);
                    i = left_bracket[left_bracket.size() - 1] - 1;
                    left_bracket.pop_back();
                }
            }

            std::string equ = equation[0];
            equation.clear();
            equation = matrix_split(equ, "()");

            expr_chain.reserve(equation.size());

            int32_t need_k = 0;
            stp_data current_i = 0;

            for (stp_data i = 0; i < equation.size(); i++)
            {
                if (equation[i] == "&")
                {
                    expr_node new_node(NodeType_Gate, "Mc", GateType_Mc, expr_chain[i - 1].Get_Var_k(), 0);
                    expr_chain.push_back(new_node);
                    need_k = 1;
                }
                else if (equation[i] == "|")
                {
                    expr_node new_node(NodeType_Gate, "Md", GateType_Md, expr_chain[i - 1].Get_Var_k(), 0);
                    expr_chain.push_back(new_node);
                    need_k = 1;
                }
                else if (equation[i] == "~")
                {
                    expr_node new_node(NodeType_Gate, "Mn", GateType_Mn, 0, 0);
                    expr_chain.push_back(new_node);
                    need_k = 1;
                }
                else if (equation[i] == "" || equation[i] == " ")
                {
                    continue;
                }
                // variable
                else 
                {
                    auto it = nameToId.find(equation[i]);
                    // existed
                    if (it != nameToId.end())
                    {
                        id index = it->second;
                        stp_data k = this->input_k[index];
                        expr_node new_node(NodeType_Variable, equation[i], GateType_NONE, index, k);
                        expr_chain.push_back(new_node);

                        if(need_k)
                        {
                            expr_chain[i - 1].Set_dim2(k);
                            need_k = 0;
                        }
                    }
                    else
                    {
                        //error
                        std::cout << "error: variable " << equation[i] << " not found" << std::endl;
                    }
                }
            }

            // print expr_chain
            if (verbose)
            {
                for (const auto &mat : expr_chain)
                {
                    std::cout << mat.Get_Node_Name() << " ";
                }
                std::cout << std::endl;
            }
        }

        void normalize_expr_chain(void)
        {
            stp_data pi_num = 0;
            expr_chain = move_vars_to_rightside(expr_chain, pi_num);
            expr_chain = sort_variables(expr_chain, pi_num);
        }

        void from_expr_to_matrix(void)
        {
            int32_t init_flag = 0;
            int32_t I_flag = 0;

            for (size_t i = 1; i < expr_chain.size(); i++)
            {
                switch (expr_chain[i].Get_GateType())
                {
                    case GateType_W:
                    {
                        if (init_flag)
                        {
                            if (I_flag)
                            {
                                matrix temp = generate_swap_matrix(expr_chain[i].Get_dim1(), expr_chain[i].Get_dim2());
                                result_matrix = Matrix_semi_tensor_product(result_matrix, In_KR_Matrix(expr_chain[i - 1].Get_dim1(), temp));
                            }
                            else
                            {
                                result_matrix = Matrix_semi_tensor_product(result_matrix, generate_swap_matrix(expr_chain[i].Get_dim1(), expr_chain[i].Get_dim2()));
                            }
                        }
                        else
                        {
                            if (I_flag)
                            {
                                matrix temp = generate_swap_matrix(expr_chain[i].Get_dim1(), expr_chain[i].Get_dim2());
                                result_matrix = Matrix_semi_tensor_product(result_matrix, In_KR_Matrix(expr_chain[i - 1].Get_dim1(), temp));
                            }
                            else
                            {
                                result_matrix = generate_swap_matrix(expr_chain[i].Get_dim1(), expr_chain[i].Get_dim2());
                            }
                            init_flag = 1;
                        }
                        break;
                    }
                    case GateType_I:
                    {
                        I_flag = 1;
                        break;
                    }
                    case GateType_Mr:
                    {
                        if (init_flag)
                        {
                            if (I_flag)
                            {
                                matrix temp = generate_Mr_matrix(expr_chain[i].Get_dim1());
                                result_matrix = Matrix_semi_tensor_product(result_matrix, In_KR_Matrix(expr_chain[i - 1].Get_dim1(), temp));
                            }
                            else
                            {
                                result_matrix = Matrix_semi_tensor_product(result_matrix, generate_Mr_matrix(expr_chain[i].Get_dim1()));
                            }
                        }
                        else
                        {
                            if (I_flag)
                            {
                                matrix temp = generate_Mr_matrix(expr_chain[i].Get_dim1());
                                result_matrix = Matrix_semi_tensor_product(result_matrix, In_KR_Matrix(expr_chain[i - 1].Get_dim1(), temp));
                            }
                            else
                            {
                                result_matrix = generate_Mr_matrix(expr_chain[i].Get_dim1());
                            }
                            init_flag = 1;
                        }
                        break;
                    }
                    case GateType_Mc:
                    {
                        if (init_flag)
                        {
                            if (I_flag)
                            {
                                matrix temp = generate_Mc_matrix(expr_chain[i].Get_Var_k());
                                result_matrix = Matrix_semi_tensor_product(result_matrix, In_KR_Matrix(expr_chain[i - 1].Get_dim1(), temp));
                            }
                            else
                            {
                                result_matrix = Matrix_semi_tensor_product(result_matrix, generate_Mr_matrix(expr_chain[i].Get_dim1()));
                            }
                        }
                        else
                        {
                            if (I_flag)
                            {
                                matrix temp = generate_Mc_matrix(expr_chain[i].Get_Var_k());
                                result_matrix = Matrix_semi_tensor_product(result_matrix, In_KR_Matrix(expr_chain[i - 1].Get_dim1(), temp));
                            }
                            else
                            {
                                result_matrix = generate_Mc_matrix(expr_chain[i].Get_Var_k());
                            }
                            init_flag = 1;
                        }
                        break;
                    }
                    case GateType_Md:
                    {
                        if (init_flag)
                        {
                            if (I_flag)
                            {
                                matrix temp = generate_Md_matrix(expr_chain[i].Get_Var_k());
                                result_matrix = Matrix_semi_tensor_product(result_matrix, In_KR_Matrix(expr_chain[i - 1].Get_dim1(), temp));
                            }
                            else
                            {
                                result_matrix = Matrix_semi_tensor_product(result_matrix, generate_Md_matrix(expr_chain[i].Get_dim1()));
                            }
                        }
                        else
                        {
                            if (I_flag)
                            {
                                matrix temp = generate_Md_matrix(expr_chain[i].Get_Var_k());
                                result_matrix = Matrix_semi_tensor_product(result_matrix, In_KR_Matrix(expr_chain[i - 1].Get_dim1(), temp));
                            }
                            else
                            {
                                result_matrix = generate_Md_matrix(expr_chain[i].Get_Var_k());
                            }
                            init_flag = 1;
                        }
                        break;
                    }
                    case GateType_Mn:
                    {
                        if (init_flag)
                        {
                            if (I_flag)
                            {
                                matrix temp = generate_Mn_matrix(expr_chain[i].Get_Var_k());
                                result_matrix = Matrix_semi_tensor_product(result_matrix, In_KR_Matrix(expr_chain[i - 1].Get_dim1(), temp));
                            }
                            else
                            {
                                result_matrix = Matrix_semi_tensor_product(result_matrix, generate_Mn_matrix(expr_chain[i].Get_Var_k()));
                            }
                        }
                        else
                        {
                            if (I_flag)
                            {
                                matrix temp = generate_Mn_matrix(expr_chain[i].Get_dim1());
                                result_matrix = Matrix_semi_tensor_product(result_matrix, In_KR_Matrix(expr_chain[i - 1].Get_dim1(), temp));
                            }
                            else
                            {
                                result_matrix = generate_Mn_matrix(expr_chain[i].Get_Var_k());
                            }
                            init_flag = 1;
                        }
                        break;
                    }
                    default:
                    {
                        std::cout << "error: gate type not found" << std::endl;
                        break;
                    }
                }
            }
        }

        // exchange variables 
        void exchange_vars_mixed(void)
        {
            int64_t new_num = 0;
            int32_t remain = 0;

            std::vector<stp_data> vec = Matrix_to_tt(result_matrix);
 
            std::vector<stp_data> out_vec(vec.size());

            for(size_t i = 0; i < out_vec.size(); i++)
            {
                new_num = 0;
                remain = i;

                for (size_t j = 0; j < old_pi_list.size(); j++)
                {
                    new_num += (remain / old_pi_list[j].Get_Var_k()) * pow(old_pi_list[j].Get_Var_k(), k - 1 - find_id(old_pi_list[j].Get_Var_id()));
                    remain = remain % old_pi_list[j].Get_Var_k();
                }

                out_vec[i] = vec[new_num];
            }


            //print tt
            for(size_t i = out_vec.size() - 1; i >= 0; i--)
            {
                std::cout << "tt:" << std::endl;
                if(out_vec[i] == 0 || out_vec[i] == 1)
                {
                    std::cout << out_vec[i] <<" "<< std::endl;
                }
                else
                {
                    std::cout << out_vec[i] << "/" <<  k-1 << " " <<std::endl;
                }
            }
        }



        //move variables to right side
        std::vector<expr_node> move_vars_to_rightside(const std::vector<expr_node> &expr_chain, stp_data &pi_num)
        {
            std::vector<expr_node> new_chain;
            std::vector<expr_node> pi_chain;
            stp_data I_dim = 1;

            //move variables to right side
            for (size_t i = 0; i < expr_chain.size(); i++)
            {
                //gate
                if (expr_chain[i].Get_NodeType() == NodeType_Gate)
                {
                    if(I_dim != 1)
                    {
                        expr_node new_node(NodeType_Gate, "I", GateType_I, I_dim);
                        new_chain.push_back(new_node);
                        I_dim = 1;
                    }
                    new_chain.push_back(expr_chain[i]);
                }
                //variable
                else
                {
                    pi_chain.push_back(expr_chain[i]);
                    I_dim *= expr_chain[i].Get_dim1();
                }
            }
            pi_num = pi_chain.size();

            new_chain.insert(new_chain.end(), pi_chain.begin(), pi_chain.end());

            return new_chain;
        }

        // sort variables
        std::vector<expr_node> sort_variables( const std::vector<expr_node> &expr_chain, const stp_data &pi_num)
        {
            std::vector<expr_node> pi_chain;
            std::vector<expr_node> new_chain;
            stp_data W_dim = 1;
            stp_data I_dim = 1;

            //copy variables to pi_chain
            for (size_t i = expr_chain.size() - pi_num; i < expr_chain.size(); i++)
            {
                pi_chain.push_back(expr_chain[i]);
            }

            for (size_t i = pi_chain.size() - 1; i >= 0; i--)
            {
                for (size_t j = i - 1; j >= 0; j--)
                {
                    // equal id
                    if (pi_chain[i].Get_Var_id() == pi_chain[j].Get_Var_id())
                    {
                        pi_chain.insert(pi_chain.begin() + i, pi_chain[i]);
                        pi_chain.erase(pi_chain.begin() + j);
                        I_dim = get_I_dim_before_var(pi_chain, j);
                        
                        expr_node new_node(NodeType_Gate, "I", GateType_I, I_dim);
                        new_chain.push_back(new_node);
                        expr_node new_node1(NodeType_Gate, "W", GateType_W, W_dim, pi_chain[i].Get_Var_k());
                        new_chain.push_back(new_node1);
                        expr_node new_node2(NodeType_Gate, "I", GateType_I, I_dim * W_dim);
                        new_chain.push_back(new_node2);
                        expr_node new_node3(NodeType_Gate, "Mr", GateType_Mr, pi_chain[i].Get_Var_k());
                        new_chain.push_back(new_node3);
                    }
                    else
                    {
                        W_dim *= pi_chain[j].Get_Var_k();
                    }
                }
                W_dim = 1;
            }

            new_pi_list.resize(input_names.size());
            old_pi_list.resize(input_names.size());

            new_pi_list[0] = pi_chain[0];
            id last_id = new_pi_list[0].Get_Var_id();
            stp_data last_len = 1;
            old_pi_list[last_id] = pi_chain[0];

            for (size_t i = 0; i < pi_chain.size(); i++)
            {
                if ( last_id != pi_chain[i].Get_Var_id() )
                {
                    new_pi_list[last_len] = pi_chain[i];
                    old_pi_list[pi_chain[i].Get_Var_id()] = pi_chain[i];
                    last_id++;
                    last_len++;
                }
            }

            return new_chain;
        }

        stp_data get_I_dim_before_var(const std::vector<expr_node> &pi_chain, const stp_data &var_i)
        {
            stp_data I_dim = 1;

            for (stp_data i = 0; i < var_i; i++)
            {
                I_dim *= pi_chain[i].Get_Var_k();
            }
            return I_dim;
        }


    

        std::vector<stp_data> Matrix_to_tt(matrix &A)
        {
            stp_data size = A.cols();
            std::vector<stp_data> Result(size);
            // operate on each column
            for (stp_data i = 0; i < A.cols(); i++)
            {
                for (stp_data j = 0; j < A.rows(); j++)
                {
                    // non-zero element
                    if (A(j, i))
                    {
                        Result[i] = k - j - 1;
                        break;
                    }
                }
            }
            return Result;
        }


        stp_data find_id(id index)
        {
            for (size_t i = 0; i < new_pi_list.size(); i++)
            {
                if (new_pi_list[i].Get_Var_id() == index)
                return i;
            }
        }




    private:


        std::vector<std::string> matrix_split(const std::string &input,
                                              const std::string &pred)
        {
            std::vector<std::string> result;
            std::string temp = "";

            unsigned count1 = input.size();
            unsigned count2 = pred.size();
            unsigned j;

            for (size_t i = 0; i < count1; i++)
            {
                for (j = 0; j < count2; j++)
                {
                    if (input[i] == pred[j])
                    {
                        break;
                    }
                }

                if (j == count2)
                {
                    temp += input[i];
                }
                else
                {
                    if (!temp.empty())
                    {
                        result.push_back(temp);
                        temp.clear();
                    }
                }
            }

            if (!temp.empty())
            {
                result.push_back(temp);
                temp.clear();
            }

            return result;
        }





        const std::vector<std::string>& input_names;
        std::unordered_map<std::string, id> nameToId;
        std::vector<expr_node> old_pi_list;
        std::vector<expr_node> new_pi_list;

        std::vector<expr_node> expr_chain;
        bool verbose;
        const std::string& expression;
        const std::vector<stp_data> &input_k;
        matrix result_matrix;
    };

} // namespace stp





