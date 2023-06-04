#include <iostream>
#include <set>
#include <tuple>
#include <chrono>

#define ME "Memory error"

using namespace std;
using namespace std::chrono;

typedef tuple <int, int> ids;
class dynamic_time_warping
{
protected:
    double* c;
    double* q;
    int len_c;
    int len_q;
    double** matrix_transforms;
    double* way;
    int** way_coords;
    int len_way;
public:
    dynamic_time_warping(const double* c, const double* q, const int len_c, const int len_q)
    {
        dynamic_time_warping::c = new double[len_c];
        if (dynamic_time_warping::c == nullptr) throw 'e';
        for (int i = 0; i < len_c; i++)
        {
            dynamic_time_warping::c[i] = c[i];
        }
        dynamic_time_warping::q = new double[len_q];
        if (dynamic_time_warping::q == nullptr) throw 'e';
        for (int i = 0; i < len_q; i++)
        {
            dynamic_time_warping::q[i] = q[i];
        }
        dynamic_time_warping::len_c = len_c;
        dynamic_time_warping::len_q = len_q;
        dynamic_time_warping::matrix_transforms = new double* [len_c];
        dynamic_time_warping::way = new double[len_c + len_q];
        dynamic_time_warping::len_way = len_c + len_q;
        if (dynamic_time_warping::way == nullptr) throw 'e';
        if (dynamic_time_warping::matrix_transforms == nullptr) throw 'e';
        for (int i = 0; i < len_c; i++)
        {
            dynamic_time_warping::matrix_transforms[i] = new double[len_q];
            if (dynamic_time_warping::matrix_transforms[i] == nullptr) throw 'e';
        }
        dynamic_time_warping::way_coords = new int* [len_c + len_q];
        if (dynamic_time_warping::way_coords == nullptr) throw 'e';
        for (int i = 0; i < len_c + len_q; i++)
        {
            dynamic_time_warping::way_coords[i] = new int[2];
            if (dynamic_time_warping::way_coords[i] == nullptr) throw 'e';
        }

    }


    dynamic_time_warping(const dynamic_time_warping& C)
    {
        dynamic_time_warping::c = new double[C.len_c];
        if (dynamic_time_warping::c == nullptr) throw 'e';
        for (int i = 0; i < C.len_c; i++)
        {
            dynamic_time_warping::c[i] = C.c[i];
        }
        dynamic_time_warping::q = new double[C.len_q];
        if (dynamic_time_warping::q == nullptr) throw 'e';
        for (int i = 0; i < C.len_q; i++)
        {
            dynamic_time_warping::q[i] = C.q[i];
        }
        dynamic_time_warping::len_c = C.len_c;
        dynamic_time_warping::len_q = C.len_q;
        dynamic_time_warping::matrix_transforms = new double* [C.len_c];
        dynamic_time_warping::way = new double[C.len_c + C.len_q];
        dynamic_time_warping::len_way = C.len_way;
        if (dynamic_time_warping::way == nullptr) throw 'e';
        if (dynamic_time_warping::matrix_transforms == nullptr) throw 'e';
        for (int i = 0; i < C.len_c; i++)
        {
            dynamic_time_warping::matrix_transforms[i] = new double[C.len_q];
            if (dynamic_time_warping::matrix_transforms[i] == nullptr) throw 'e';
        }
        for (int i = 0; i < C.len_c; i++)
        {
            for (int j = 0; j < C.len_q; j++)
            {
                dynamic_time_warping::matrix_transforms[i][j] = C.matrix_transforms[i][j];
            }
        }
        dynamic_time_warping::way_coords = new int* [C.len_way];
        if (dynamic_time_warping::way_coords == nullptr) throw 'e';
        for (int i = 0; i < C.len_way; i++)
        {
            dynamic_time_warping::way_coords[i] = new int[2];
            if (dynamic_time_warping::way_coords[i] == nullptr) throw 'e';
        }
        for (int i = 0; i < C.len_c + C.len_q; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                dynamic_time_warping::way_coords[i][j] = C.way_coords[i][j];
            }
        }

    }

    ~dynamic_time_warping()
    {
        if (dynamic_time_warping::c != nullptr) dynamic_time_warping::c = nullptr;
        if (dynamic_time_warping::q != nullptr) dynamic_time_warping::q = nullptr;
        if (dynamic_time_warping::matrix_transforms != nullptr) dynamic_time_warping::matrix_transforms = nullptr;
        if (dynamic_time_warping::way != nullptr) dynamic_time_warping::way = nullptr;
    }

    double** get_matrix_transforms()
    {
        return dynamic_time_warping::matrix_transforms;
    }

    int** get_way()
    {
        return dynamic_time_warping::way_coords;
    }
    
    int get_len_way()
    {
        return dynamic_time_warping::len_way;
    }

    double get_last()
    {
        return dynamic_time_warping::matrix_transforms[len_c - 1][len_q - 1];
    }

    void find_matrix_transforms_default()
    {
        for (int i = 0; i < dynamic_time_warping::len_c; i++)
        {
            for (int j = 0; j < dynamic_time_warping::len_q; j++)
            {
                double d = abs(dynamic_time_warping::c[i] - dynamic_time_warping::q[j]);
                if (i == 0 && j == 0) dynamic_time_warping::matrix_transforms[i][j] = d;
                else if (i == 0 && j != 0) dynamic_time_warping::matrix_transforms[i][j] = d + dynamic_time_warping::matrix_transforms[i][j - 1];
                else if (i != 0 && j == 0) dynamic_time_warping::matrix_transforms[i][j] = d + dynamic_time_warping::matrix_transforms[i - 1][j];
                else
                {
                    dynamic_time_warping::matrix_transforms[i][j] = d +
                        dynamic_time_warping::matrix_transforms[i - 1][j] * (dynamic_time_warping::matrix_transforms[i - 1][j] < dynamic_time_warping::matrix_transforms[i][j - 1] && dynamic_time_warping::matrix_transforms[i - 1][j] < dynamic_time_warping::matrix_transforms[i - 1][j - 1]) +
                        dynamic_time_warping::matrix_transforms[i][j - 1] * (dynamic_time_warping::matrix_transforms[i][j - 1] <= dynamic_time_warping::matrix_transforms[i - 1][j] && dynamic_time_warping::matrix_transforms[i][j - 1] < dynamic_time_warping::matrix_transforms[i - 1][j - 1]) +
                        dynamic_time_warping::matrix_transforms[i - 1][j - 1] * (dynamic_time_warping::matrix_transforms[i - 1][j - 1] <= dynamic_time_warping::matrix_transforms[i - 1][j] && dynamic_time_warping::matrix_transforms[i - 1][j - 1] <= dynamic_time_warping::matrix_transforms[i][j - 1]);
                }
            } //TODO Потестить на разных данных
        }
    }

    void make_way_default()
    {
        int way_l = 0;
        int i = dynamic_time_warping::len_c - 1;
        int j = dynamic_time_warping::len_q - 1;
        while (i != 0 || j != 0)
        {
            double d = abs(dynamic_time_warping::c[i] - dynamic_time_warping::q[j]);
            dynamic_time_warping::way[way_l] = dynamic_time_warping::matrix_transforms[i][j];
            int* a = new int[2];
            a[0] = i;
            a[1] = j;
            dynamic_time_warping::way_coords[way_l] = a;
            if (dynamic_time_warping::matrix_transforms[i - 1][j - 1] + d == dynamic_time_warping::matrix_transforms[i][j])
            {
                i--;
                j--;
            }
            else if (dynamic_time_warping::matrix_transforms[i - 1][j] + d == dynamic_time_warping::matrix_transforms[i][j])
            {
                i--;
            }
            else if (dynamic_time_warping::matrix_transforms[i][j - 1] + d == dynamic_time_warping::matrix_transforms[i][j])
            {
                j--;
            }
            way_l++;
        }
        dynamic_time_warping::way[way_l] = dynamic_time_warping::matrix_transforms[0][0];
        dynamic_time_warping::len_way = way_l + 1;
        int a[2] = { 0, 0 };
        dynamic_time_warping::way_coords[way_l] = a;
    }

    void print_matrix_transforms()
    {
        cout << dynamic_time_warping::len_c << " " << dynamic_time_warping::len_q << "\n";
        for (int i = 0; i < dynamic_time_warping::len_c; i++)
        {
            for (int j = 0; j < dynamic_time_warping::len_q; j++)
            {
                cout << dynamic_time_warping::matrix_transforms[i][j] << " ";
            }
            cout << "\n";
        }
    }

    void print_way()
    {
        cout << dynamic_time_warping::len_way << "\n";
        for (int i = dynamic_time_warping::len_way - 1; i >= 0; i--) cout << "(" << dynamic_time_warping::way_coords[i][0] << ", " << dynamic_time_warping::way_coords[i][1] << ") ";
        cout << "\n";
    }

    void print()
    {
        dynamic_time_warping::print_matrix_transforms();
        dynamic_time_warping::print_way();
        cout << dynamic_time_warping::get_last();
        cout << "\n";
    }
};

class fast_dynamic_time_warping : public dynamic_time_warping
{
protected:
    int radius;
    ids* window;
    int len_window;
public:
    fast_dynamic_time_warping(const double* c, const double* q, const int len_c, const int len_q, const int radius) : dynamic_time_warping(c, q, len_c, len_q)
    {
        fast_dynamic_time_warping::radius = radius;
        fast_dynamic_time_warping::window = new ids[len_q * len_c];
        fast_dynamic_time_warping::len_window = 0;
    }

    fast_dynamic_time_warping(const fast_dynamic_time_warping& C) : dynamic_time_warping(C.c, C.q, C.len_c, C.len_q)
    {
        fast_dynamic_time_warping::radius = C.radius;
        fast_dynamic_time_warping::window = new ids[C.len_q * C.len_c];
        fast_dynamic_time_warping::len_window = 0;
    }

    ~fast_dynamic_time_warping() {}
    void fast_dtw()
    {
       fast_dtw_inner(fast_dynamic_time_warping::c, fast_dynamic_time_warping::q, fast_dynamic_time_warping::len_c, fast_dynamic_time_warping::len_q, fast_dynamic_time_warping::radius);
    }
    
    void fast_dtw_inner(const double* c, const double* q, const int len_c, const int len_q, const int radius)
    {
        if (len_c < radius + 2 || len_q < radius + 2)
        {
            if (window == nullptr) throw 'e';
            len_window = 0;
            for (int i = 0; i < len_c; i++)
            {
                for (int j = 0; j < len_q; j++)
                {
                    window[len_window] = ids(i, j);
                    len_window++;
                    
                }
            }
            local_dtw(c, q, len_c, len_q);
            return;
        }
        double* c_shrinked = reduce_by_half(c, len_c);
        if (c_shrinked == nullptr) throw 'e';
        double* q_shrinked = reduce_by_half(q, len_q);
        if (q_shrinked == nullptr) throw 'e';
        fast_dtw_inner(c_shrinked, q_shrinked, len_c / 2, len_q / 2, radius);
        expand_window(len_c, len_q, radius);
        local_dtw(c, q, len_c, len_q);
        return;
    }

    double* reduce_by_half(const double* x, const int len)
    {
        double* local_x = new double[len/2];
        if (local_x == nullptr) throw 'e';
        for (int i = 0; i < (len / 2); i++)
        {
            local_x[i] = (x[i * 2] + x[i * 2 + 1]) / 2;
        }
        return local_x;
    }

    void expand_window(const int len_c, const int len_q, const int radius) 
    {
        std::set<ids> local_path;
        for (int i = 0; i < len_way; i++)
        {
            local_path.insert(ids(way_coords[i][0], way_coords[i][1]));
        }
        for (int i = 0; i < len_way; i++)
        {
            for (int a = -radius; a < radius + 1; a++)
            {
                for (int b = -radius; b < radius + 1; b++)
                {
                    local_path.insert(ids(a, b));
                }
            }
        }
        std::set<ids> local_window; // Мб считать длину?
        for (auto itr = local_path.begin(); itr != local_path.end(); itr++)
        {
            ids b = static_cast<ids>(*itr);
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    local_window.insert(ids(get<0>(b) * 2 + i, get<1>(b) * 2 + j));
                }
            }
        }
        len_window = 0;
        int start_j = 0;
        for (int i = 0; i < len_c; i++)
        {
            int new_start_j = NULL;
            for (int j = start_j; j < len_q; j++)
            {
                ids curr = ids(i, j);
                if (local_window.count(curr) != 0)
                {
                    window[len_window] = curr;
                    len_window++;
                    if (!new_start_j)
                    {
                        new_start_j = j;
                    }
                }
                else if (new_start_j)
                {
                    break;
                }
            }
            start_j = new_start_j;
        }
    }
    void local_dtw(const double* c, const double* q, const int dlen_c, const int dlen_q)
    {
        double** loc_matrix = new double*[len_c];
        if (loc_matrix == nullptr) throw 'e';
        for (int i = 0; i < dlen_c; i++)
        {
            loc_matrix[i] = new double[len_q];
            if (loc_matrix[i] == nullptr) throw 'e';
            for (int j = 0; j < dlen_q; j++)
            {
                loc_matrix[i][j] = INFINITY;
            }
        }
        for (int i = 0; i < len_window; i++)
        {
            int a = get<0>(window[i]);
            int b = get<1>(window[i]);
            loc_matrix[a][b] = abs(q[a] - c[b]);
            if (a != 0 && b != 0)
            {
                loc_matrix[a][b] += min(min(loc_matrix[a][b - 1], loc_matrix[a - 1][b]), loc_matrix[a - 1][b - 1]);
            }
            else if (a != 0 && b == 0)
            {
                loc_matrix[a][b] += loc_matrix[a - 1][b];
            }
            else if (a == 0 && b != 0)
            {
                loc_matrix[a][b] += loc_matrix[a][b - 1];
            }
        }
        int** loc_way = new int*[dlen_c + dlen_q];
        if (loc_way == nullptr) throw 'e';
        int len_loc_way = 0;
        int i = dlen_c - 1;
        int j = dlen_q - 1;
        while (i != 0 || j != 0)
        {
            loc_way[len_loc_way] = new int[2];
            if (loc_way[len_loc_way] == nullptr) throw 'e';
            loc_way[len_loc_way][0] = j;
            loc_way[len_loc_way][1] = i;
            len_loc_way++;
            if (i != 0 && j != 0)
            {
                if (min(min(loc_matrix[i][j - 1], loc_matrix[i - 1][j]), loc_matrix[i - 1][j - 1]) == loc_matrix[i - 1][j - 1])
                {
                    i--;
                    j--;
                }
                else if (min(min(loc_matrix[i][j - 1], loc_matrix[i - 1][j]), loc_matrix[i - 1][j - 1]) == loc_matrix[i][j - 1])
                {
                    j--;
                }
                else
                {
                    i--;
                }
            }
            else if (i == 0 && j != 0)
            {
                j--;
            }
            else if (i != 0 && j == 0)
            {
                i--;
            }
        }
        
        len_way = 0;
        for (int i = 0; i < len_loc_way; i++)
        {
            way_coords[len_way] = new int[2];
            if (way_coords[len_way] == nullptr) throw 'e';
            way_coords[len_way][0] = loc_way[i][0];
            way_coords[len_way][1] = loc_way[i][1];
            len_way++;
        }
        way_coords[len_way] = new int[2];
        if (way_coords[len_way] == nullptr) throw 'e';
        way_coords[len_way][0] = 0;
        way_coords[len_way][1] = 0;
        len_way++;
    }
};



class Sparce_DTW : public  dynamic_time_warping
{
protected:
    double res;
    double* quant_q;
    double* quant_c;
public:

    // Конструктор, пока только заполненый, потом допишу пустой
    Sparce_DTW(double* c_seq, double* q_seq, const int c_l, const int q_l, double r) : dynamic_time_warping(c_seq, q_seq, c_l, q_l) {
        if (!(quant_q = new double[q_l])) { throw ME; }
        if (!(quant_c = new double[c_l])) { throw ME; }
        res = r;

        quant_q = quantanize(q, len_q);
        quant_c = quantanize(c, len_c);
    }

    ~Sparce_DTW()
    {
        if (Sparce_DTW::c != nullptr) Sparce_DTW::c = nullptr;
        if (Sparce_DTW::q != nullptr) Sparce_DTW::q = nullptr;
        if (Sparce_DTW::quant_c != nullptr) Sparce_DTW::c = nullptr;
        if (Sparce_DTW::quant_q != nullptr) Sparce_DTW::q = nullptr;
        if (Sparce_DTW::matrix_transforms != nullptr) Sparce_DTW::matrix_transforms = nullptr;
        if (Sparce_DTW::way != nullptr) Sparce_DTW::way = nullptr;
    }

    void ret_q_quant() const
    {
        for (int i = 0; i < len_q; i++)
        {
            cout << quant_q[i] << " ";
        }
        cout << "\n";
    }

    void ret_c_quant() const
    {
        for (int i = 0; i < len_c; i++)
        {
            cout << quant_c[i] << " ";
        }
        cout << "\n";
    }

    double* quantanize(double* seq, int seq_l)
    {
        double* new_quant_seq;
        if (!(new_quant_seq = new double[seq_l])) { throw ME; }
        double mn = 100000;
        double mx = -100000;
        for (int i = 0; i < seq_l; i++)
        {
            if (seq[i] > mx) { mx = seq[i]; }
            if (seq[i] < mn) { mn = seq[i]; }
        }
        for (int i = 0; i < seq_l; i++)
        {
            new_quant_seq[i] = (seq[i] - mn) / (mx - mn);
        }

        return new_quant_seq;
    }

    double euc_dist(double x, double y)
    {
        return abs(x - y);
    }

    void matrix_transforms_to_zeros()
    {
        for (int i = 0; i < len_c; i++)
        {
            for (int j = 0; j < len_q; j++)
            {
                matrix_transforms[i][j] = 0;
            }
        }

    }

    void populate_warp()
    {
        matrix_transforms_to_zeros();
        double lower_bound = 0;
        double upper_bound = res;
        while (0 <= lower_bound && (1 - (res / 2)) >= lower_bound)
        {
            int* idxC;
            int* idxQ;

            if (!(idxC = new int[len_c])) { throw ME; }
            if (!(idxQ = new int[len_q])) { throw ME; }

            for (int i = 0; i < len_c; i++)
            {
                if (lower_bound <= quant_c[i] && quant_c[i] <= upper_bound) { idxC[i] = 1; }
            }

            for (int i = 0; i < len_q; i++)
            {
                if (lower_bound <= quant_q[i] && quant_q[i] <= upper_bound) { idxQ[i] = 1; }
            }

            lower_bound = lower_bound + res / 2;
            upper_bound = lower_bound + res;

            for (int i = 0; i < len_q; i++)
            {
                for (int j = 0; j < len_c; j++)
                {
                    double euc_d;
                    if (idxQ[i] == 1 && idxC[j] == 1)
                    {
                        euc_d = euc_dist(quant_c[j], quant_q[i]);
                        if (euc_d == 0)
                        {
                            matrix_transforms[j][i] = -1;
                        }
                        else
                        {
                            if (matrix_transforms[j][i] != -1)
                            {
                                matrix_transforms[j][i] = matrix_transforms[j][i] + euc_d;
                            }
                            else
                            {
                                matrix_transforms[j][i] = euc_d;
                            }
                        }
                    }
                }
            }

        }
    }

    int** lower_neighbors(int x, int y)
    {
        int** ln;
        if (!(ln = new int* [3])) { throw ME; }
        for (int i = 0; i < 3; i++)
        {
            if (!(ln[i] = new int[2])) { throw ME; }
        }
        if (x == 0 && y == 0)
        {
            return NULL;
        }
        else
        {
            if (x - 1 >= 0 && y >= 0)
            {
                ln[0][0] = x - 1;
                ln[0][1] = y;
            }
            else
            {
                ln[0] = NULL;
            }

            if (x - 1 >= 0 && y - 1 >= 0)
            {
                ln[1][0] = x - 1;
                ln[1][1] = y - 1;
            }
            else { ln[1] = NULL; }

            if (x >= 0 && y - 1 >= 0)
            {
                ln[2][0] = x;
                ln[2][1] = y - 1;
            }
            else { ln[2] = NULL; }

            return ln;
        }
    }


    int** upper_neighbors(int x, int y, int max_x, int max_y)
    {
        int** un;
        if (!(un = new int* [3])) { throw ME; }
        for (int i = 0; i < 3; i++)
        {
            if (!(un[i] = new int[2])) { throw ME; }
        }
        if (x == max_x && y == max_y)
        {
            return NULL;
        }
        else
        {
            if (x + 1 <= max_x && y <= max_y)
            {
                un[0][0] = x + 1;
                un[0][1] = y;
            }
            else { un[0] = NULL; }

            if (x + 1 <= max_x && y + 1 <= max_y)
            {
                un[1][0] = x + 1;
                un[1][1] = y + 1;
            }
            else { un[1] = NULL; }

            if (x <= max_x && y + 1 <= max_y)
            {
                un[2][0] = x;
                un[2][1] = y + 1;
            }
            else { un[2] = NULL; }

            return un;
        }
    }

    void unlock_upper_neighbors(int** neighbors)
    {
        for (int i = 0; i < 3; i++)
        {
            int x = neighbors[i][0];
            int y = neighbors[i][1];
            if (matrix_transforms[x][y] == 0)
            {
                matrix_transforms[x][y] = euc_dist(quant_c[x], quant_q[y]);
            }
        }
    }


    void calculate_warp_costs()
    {
        for (int i = 0; i < len_c; i++)
        {
            for (int j = 0; j < len_q; j++)
            {
                if (matrix_transforms[i][j])
                {
                    int** lower_n;
                    if (!(lower_n = new int* [3])) { throw ME; }
                    lower_n = lower_neighbors(i, j);
                    double min_cost = 100000;
                    if (lower_n != NULL)
                    {
                        for (int x = 0; x < 3; x++)
                        {
                            if (lower_n[x] != NULL)
                            {
                                if (matrix_transforms[lower_n[x][0]][lower_n[x][1]] == 0) { lower_n[x] = NULL; }
                            }
                        }
                        for (int x = 0; x < 3; x++)
                        {
                            if (lower_n[x] != NULL)
                            {
                                if (matrix_transforms[lower_n[x][0]][lower_n[x][1]] < min_cost) { min_cost = matrix_transforms[lower_n[x][0]][lower_n[x][1]]; }
                                if (min_cost == -1) { min_cost = 0; }
                            }
                        }
                    }
                    else { min_cost = matrix_transforms[i][j]; }

                    if (matrix_transforms[i][j] > -1) { matrix_transforms[i][j] = matrix_transforms[i][j] + min_cost; }
                    else if (matrix_transforms[i][j] == -1 && min_cost > 0) { matrix_transforms[i][j] = min_cost; }

                    int** upper_n;
                    if (!(upper_n = new int* [3])) { throw ME; }
                    upper_n = upper_neighbors(i, j, len_c - 1, len_q - 1);
                    if (upper_n != NULL)
                    {
                        for (int x = 0; x < 3; x++)
                        {
                            if (upper_n[x] != NULL)
                            {
                                if (matrix_transforms[upper_n[x][0]][upper_n[x][1]] == 0) { break; }
                                if (x == 2) { unlock_upper_neighbors(upper_n); }
                            }
                        }
                    }
                }
            }
        }

    }

    void calculate_warp_path()
    {
        int x = len_c - 1;
        int y = len_q - 1;

        way_coords[0][0] = x;
        way_coords[0][1] = y;

        int f = 1;

        while (x != 0 && y != 0)
        {
            int** lower_n;
            if (!(lower_n = new int* [3])) { throw ME; }
            lower_n = lower_neighbors(x, y);

            for (int i = 0; i < 3; i++)
            {
                if (lower_n[i] != NULL)
                {
                    if (matrix_transforms[lower_n[i][0]][lower_n[i][1]] == 0) { lower_n[i] = NULL; }
                }
            }
            if (lower_n[0] == NULL && lower_n[1] == NULL && lower_n[2] == NULL) { lower_n = NULL; }

            if (lower_n != NULL)
            {


                int* lowest;
                if (!(lowest = new int[2])) { throw ME; }


                lowest = lower_n[0];
                for (int i = 1; i < 3; i++)
                {
                    if (lowest == NULL) { lowest = lower_n[i]; }
                    else if (i < 3 && lower_n[i] != NULL)
                    {
                        if (matrix_transforms[lowest[0]][lowest[1]] > matrix_transforms[lower_n[i][0]][lower_n[i][1]]) { lowest = lower_n[i]; }
                    }
                }
                way_coords[f][0] = lowest[0];
                way_coords[f][1] = lowest[1];
                x = lowest[0];
                y = lowest[1];
                f++;
            }
        }
        way_coords[f][0] = 0;
        way_coords[f][1] = 0;
        len_way = f + 1;
        for (int i = 0; i < len_way; i++)
        {
            way[i] = matrix_transforms[way_coords[i][0]][way_coords[i][1]];
        }
    }
};


int main()
{
    double a[16] = { 1, 2, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    double b[16] = { 1, 5, 4, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    auto s = steady_clock::now();
    dynamic_time_warping w = dynamic_time_warping(a, b, 16, 16);
    w.find_matrix_transforms_default();
    w.make_way_default();
    w.print_way();
    auto f = steady_clock::now();

    cout << "Time for cDTW " << duration_cast<chrono::nanoseconds>(f - s).count() << "\n";

    s = steady_clock::now();

    fast_dynamic_time_warping fw = fast_dynamic_time_warping(a, b, 16, 16, 1);
    fw.fast_dtw();
    fw.print_way();

    f = steady_clock::now();

    cout << "Time for FAST-DTW " << duration_cast<chrono::nanoseconds>(f - s).count() << "\n";

    s = steady_clock::now();

    Sparce_DTW sw = Sparce_DTW(a, b, 16, 16, 0.5);
    sw.populate_warp();
    sw.calculate_warp_costs();
    sw.calculate_warp_path();
    sw.print_way();

    f = steady_clock::now();

    cout << "Time for SPARSE-DTW " << duration_cast<chrono::nanoseconds>(f - s).count() << "\n";
}


// histogram equalization 
// x - double, y - int will be more reliable
