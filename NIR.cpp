#include <iostream>

using namespace std;

class dynamic_time_warping
{
protected:
    int* c;
    int* q;
    int len_c;
    int len_q;
    int** matrix_transforms;
    int* way;
    int** way_coords;
    int len_way;
public:
    dynamic_time_warping(const int* c, const int* q, const int len_c, const int len_q) 
    { 
        dynamic_time_warping::c = new int[len_c];
        if (dynamic_time_warping::c == nullptr) throw 'e';
        for (int i = 0; i < len_c; i++)
        {
            dynamic_time_warping::c[i] = c[i];
        }
        dynamic_time_warping::q = new int[len_q];
        if (dynamic_time_warping::q == nullptr) throw 'e';
        for (int i = 0; i < len_q; i++)
        {
            dynamic_time_warping::q[i] = q[i];
        }
        dynamic_time_warping::len_c = len_c;
        dynamic_time_warping::len_q = len_q;
        dynamic_time_warping::matrix_transforms = new int * [len_c];
        dynamic_time_warping::way = new int [len_c + len_q];
        dynamic_time_warping::len_way = len_c + len_q;
        if (dynamic_time_warping::way == nullptr) throw 'e';
        if (dynamic_time_warping::matrix_transforms == nullptr) throw 'e';
        for (int i = 0; i < len_c; i++)
        {
            dynamic_time_warping::matrix_transforms[i] = new int[len_q];
            if (dynamic_time_warping::matrix_transforms[i] == nullptr) throw 'e';
        }
        dynamic_time_warping::way_coords = new int*[len_c + len_q];
        if (dynamic_time_warping::way_coords == nullptr) throw 'e';
        for (int i = 0; i < len_c + len_q; i++)
        {
            dynamic_time_warping::way_coords[i] = new int[2];
            if (dynamic_time_warping::way_coords[i] == nullptr) throw 'e';
        }

    }


    dynamic_time_warping(const dynamic_time_warping& C)
    {
        dynamic_time_warping::c = new int[C.len_c];
        if (dynamic_time_warping::c == nullptr) throw 'e';
        for (int i = 0; i < C.len_c; i++)
        {
            dynamic_time_warping::c[i] = C.c[i];
        }
        dynamic_time_warping::q = new int[C.len_q];
        if (dynamic_time_warping::q == nullptr) throw 'e';
        for (int i = 0; i < C.len_q; i++)
        {
            dynamic_time_warping::q[i] = C.q[i];
        }
        dynamic_time_warping::len_c = C.len_c;
        dynamic_time_warping::len_q = C.len_q;
        dynamic_time_warping::matrix_transforms = new int* [C.len_c];
        dynamic_time_warping::way = new int[C.len_c + C.len_q];
        dynamic_time_warping::len_way = C.len_way;
        if (dynamic_time_warping::way == nullptr) throw 'e';
        if (dynamic_time_warping::matrix_transforms == nullptr) throw 'e';
        for (int i = 0; i < C.len_c; i++)
        {
            dynamic_time_warping::matrix_transforms[i] = new int[C.len_q];
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

    int** get_matrix_transforms()
    {
        return dynamic_time_warping::matrix_transforms;
    }

    int* get_way()
    {
        return dynamic_time_warping::way;
    }

    int get_last()
    {
        return dynamic_time_warping::matrix_transforms[len_c - 1][len_q - 1];
    }

    void find_matrix_transforms_default()
    {
        for (int i = 0; i < dynamic_time_warping::len_c; i++)
        {
            for (int j = 0; j < dynamic_time_warping::len_q; j++)
            {
                int d = abs(dynamic_time_warping::c[i] - dynamic_time_warping::q[j]);
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
            int d = abs(dynamic_time_warping::c[i] - dynamic_time_warping::q[j]);
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
        dynamic_time_warping::len_way = way_l+1;
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
        for (int i = dynamic_time_warping::len_way - 1; i >= 0; i--) cout << "(" << dynamic_time_warping::way_coords[i][0] << ", " << dynamic_time_warping::way_coords[i][1] << ") " << dynamic_time_warping::way[i] << " ";
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


int main()
{
    int a[4] = {1, 2, 4, 1};
    int b[4] = {1, 5, 4, 2};
    dynamic_time_warping w = dynamic_time_warping(a, b, 4, 4);
    w.find_matrix_transforms_default();
    w.print_matrix_transforms();
    w.make_way_default();
    w.print_way();
}

