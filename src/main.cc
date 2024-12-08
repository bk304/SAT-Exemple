#include <bits/stdc++.h>

#include "minisat/core/Solver.h"
#include "minisat/core/SolverTypes.h"
using namespace std;

/*
As regras que definem a evolução de um estado para outro do jogo são:
  - toda célula morta com exatamente três vizinhos vivos torna-se viva;
  - toda célula viva com menos de dois vizinhos vivos morre;
  - toda célula viva com mais de três vizinhos vivos morre;
  - toda célula viva com dois ou três vizinhos vivos permanece viva.
*/

void addUnsatResult(
    int n,
    int m,
    std::map<int, vector<vector<int>>>& resultBoards
) {
    vector<vector<int>> result(n - 2, vector<int>(m - 2, 0));
    resultBoards[-1] = result;
}

void vector_to_MinisatVec(
    vector<Minisat::Lit>& v,
    Minisat::vec<Minisat::Lit>& result
) {
    for (auto& lit : v) {
        result.push(lit);
    }
}

void generateSubsets(
    const vector<Minisat::Lit>& A,
    int k,
    int index,
    vector<Minisat::Lit>& current,
    vector<vector<Minisat::Lit>>& result
) {
    if (current.size() == k) {
        result.push_back(current);
        return;
    }

    for (int i = index; i < A.size(); ++i) {
        current.push_back(A[i]);
        generateSubsets(A, k, i + 1, current, result);
        current.pop_back();
    }
}

vector<vector<Minisat::Lit>>
subsetsOfSizeK(const vector<Minisat::Lit>& A, int k) {
    vector<vector<Minisat::Lit>> result;
    vector<Minisat::Lit> current;
    generateSubsets(A, k, 0, current, result);
    return result;
}

void encodeRelaxedDeadCells(
    Minisat::Solver& solver,
    vector<Minisat::Lit> const& cell_lits,
    vector<Minisat::Lit> const& relaxation_lits,
    int k
) {
    for (int i = 0; i < cell_lits.size(); ++i) {
        Minisat::vec<Minisat::Lit> clause;
        clause.push(~(cell_lits[i]));
        if (k > 0)
            clause.push(relaxation_lits[i]);
        solver.addClause(clause);
    }
}

// Limitarei a quantidade máxima de variáveis relaxadas habilitadas usando
// sequential counter encoding, criado por Sinz.
void enforceRelaxationBound(
    Minisat::Solver& solver,
    vector<Minisat::Var> const& relaxationVars,
    int k
) {
    if (k == 0)
        return;

    int n = relaxationVars.size();

    vector<vector<Minisat::Var>> counterRegisterSet(n);
    for (int i = 0; i < n; ++i) {
        counterRegisterSet[i].resize(k);
        for (int j = 0; j < k; ++j) {
            counterRegisterSet[i][j] = solver.newVar(false, true);
        }
    }

    {
        Minisat::vec<Minisat::Lit> clause;
        clause.push(~Minisat::mkLit(relaxationVars[0]));
        clause.push(Minisat::mkLit(counterRegisterSet[0][0]));
        solver.addClause(clause);
    }

    for (int j = 1; j < k; ++j) {
        Minisat::vec<Minisat::Lit> clause;
        clause.push(~Minisat::mkLit(counterRegisterSet[0][j]));
        solver.addClause(clause);
    }

    for (int i = 1; i < n - 1; ++i) {
        {
            Minisat::vec<Minisat::Lit> clause;
            clause.push(~Minisat::mkLit(relaxationVars[i]));
            clause.push(Minisat::mkLit(counterRegisterSet[i][0]));
            solver.addClause(clause);
        }

        {
            Minisat::vec<Minisat::Lit> clause;
            clause.push(~Minisat::mkLit(counterRegisterSet[i - 1][0]));
            clause.push(Minisat::mkLit(counterRegisterSet[i][0]));
            solver.addClause(clause);
        }

        for (int j = 1; j < k; ++j) {
            {
                Minisat::vec<Minisat::Lit> clause;
                clause.push(~Minisat::mkLit(relaxationVars[i]));
                clause.push(~Minisat::mkLit(counterRegisterSet[i - 1][j - 1]));
                clause.push(Minisat::mkLit(counterRegisterSet[i][j]));
                solver.addClause(clause);
            }

            {
                Minisat::vec<Minisat::Lit> clause;
                clause.push(~Minisat::mkLit(counterRegisterSet[i - 1][j]));
                clause.push(Minisat::mkLit(counterRegisterSet[i][j]));
                solver.addClause(clause);
            }
        }

        {
            Minisat::vec<Minisat::Lit> clause;
            clause.push(~Minisat::mkLit(relaxationVars[i]));
            clause.push(~Minisat::mkLit(counterRegisterSet[i - 1][k - 1]));
            solver.addClause(clause);
        }
    }

    {
        Minisat::vec<Minisat::Lit> clause;
        clause.push(~Minisat::mkLit(relaxationVars[n - 1]));
        clause.push(~Minisat::mkLit(counterRegisterSet[n - 2][k - 1]));
        solver.addClause(clause);
    }
}

void loneliness(
    Minisat::Solver& solver,
    vector<Minisat::Lit> const& neighbors
) {
    vector<vector<Minisat::Lit>> A = subsetsOfSizeK(neighbors, 7);
    for (auto& subset : A) {
        Minisat::vec<Minisat::Lit> subset_minisatVec;
        vector_to_MinisatVec(subset, subset_minisatVec);

        solver.addClause(subset_minisatVec);
    }
}

void stagnation(
    Minisat::Solver& solver,
    vector<Minisat::Lit> const& neighbors,
    Minisat::Lit const& cell
) {
    vector<vector<Minisat::Lit>> A = subsetsOfSizeK(neighbors, 2);
    for (auto& subset : A) {
        Minisat::vec<Minisat::Lit> subset_minisatVec;
        vector_to_MinisatVec(subset, subset_minisatVec);

        Minisat::vec<Minisat::Lit> clause;
        clause.push(cell);

        for (size_t k = 0; k < subset_minisatVec.size(); ++k) {
            clause.push(~subset_minisatVec[k]);
        }

        Minisat::vec<Minisat::Lit> temp;
        std::unordered_set<int> c_values;
        for (int i = 0; i < subset_minisatVec.size(); ++i) {
            c_values.insert(subset_minisatVec[i].x);
        }
        for (int i = 0; i < neighbors.size(); ++i) {
            if (c_values.find(neighbors[i].x)
                == c_values.end()) { // Se não está em c
                temp.push(neighbors[i]);
            }
        }
        for (size_t k = 0; k < temp.size(); ++k) {
            clause.push(temp[k]);
        }

        solver.addClause(clause);
    }
}

void overcrowding(
    Minisat::Solver& solver,
    vector<Minisat::Lit> const& neighbors
) {
    vector<vector<Minisat::Lit>> A = subsetsOfSizeK(neighbors, 4);
    for (auto& subset : A) {
        Minisat::vec<Minisat::Lit> subset_minisatVec;
        vector_to_MinisatVec(subset, subset_minisatVec);

        Minisat::vec<Minisat::Lit> clause;
        for (size_t k = 0; k < subset_minisatVec.size(); ++k) {
            clause.push(~subset_minisatVec[k]);
        }
        solver.addClause(clause);
    }
}

void preservation(
    Minisat::Solver& solver,
    vector<Minisat::Lit> const& neighbors,
    Minisat::Lit const& cell
) {
    vector<vector<Minisat::Lit>> A = subsetsOfSizeK(neighbors, 2);
    for (auto& subset : A) {
        Minisat::vec<Minisat::Lit> subset_minisatVec;
        vector_to_MinisatVec(subset, subset_minisatVec);

        Minisat::vec<Minisat::Lit> clause;
        clause.push(~cell);

        for (size_t k = 0; k < subset_minisatVec.size(); ++k) {
            clause.push(~subset_minisatVec[k]);
        }

        Minisat::vec<Minisat::Lit> temp;
        std::unordered_set<int> c_values;
        for (int i = 0; i < subset_minisatVec.size(); ++i) {
            c_values.insert(subset_minisatVec[i].x);
        }
        for (int i = 0; i < neighbors.size(); ++i) {
            if (c_values.find(neighbors[i].x)
                == c_values.end()) { // Se não está em c
                temp.push(neighbors[i]);
            }
        }
        for (size_t k = 0; k < temp.size(); ++k) {
            clause.push(temp[k]);
        }

        solver.addClause(clause);
    }
}

void life(Minisat::Solver& solver, vector<Minisat::Lit> const& neighbors) {
    vector<vector<Minisat::Lit>> A = subsetsOfSizeK(neighbors, 3);
    for (auto& subset : A) {
        Minisat::vec<Minisat::Lit> subset_minisatVec;
        vector_to_MinisatVec(subset, subset_minisatVec);

        Minisat::vec<Minisat::Lit> clause;

        for (size_t k = 0; k < subset_minisatVec.size(); ++k) {
            clause.push(~subset_minisatVec[k]);
        }

        Minisat::vec<Minisat::Lit> temp;
        std::unordered_set<int> c_values;
        for (int i = 0; i < subset_minisatVec.size(); ++i) {
            c_values.insert(subset_minisatVec[i].x);
        }
        for (int i = 0; i < neighbors.size(); ++i) {
            if (c_values.find(neighbors[i].x)
                == c_values.end()) { // Se não está em c
                temp.push(neighbors[i]);
            }
        }
        for (size_t k = 0; k < temp.size(); ++k) {
            clause.push(temp[k]);
        }

        solver.addClause(clause);
    }
}

bool solve(
    vector<vector<int>>& board,
    int k,
    std::map<int, vector<vector<int>>>& resultBoards
) {
    int n = board.size();
    int m = board[0].size();
    Minisat::Solver solver;

    // criarei uma clausula para cada celula que obriga ela a estar morta caso a variavel de relaxamento seja falso.
    // a ideia é ir aumentando a quantidade de variaveis relaxadas até que haja um resultado SAT.
    vector<Minisat::Var> relaxationVars;
    vector<Minisat::Lit> cellsLits;

    // Cria uma variável do MiniSat para cada celula do tabuleiro.
    vector<vector<Minisat::Var>> previous(n, vector<Minisat::Var>(m));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            previous[i][j] =
                solver.newVar((i == 0 || i == n - 1 || j == 0 || j == m - 1));
            relaxationVars.push_back(solver.newVar(false, false));
            cellsLits.push_back(Minisat::mkLit(previous[i][j]));
        }
    }

    // Cria as clausulas que limitam a quantidade de celulas vivas usando variáveis relaxadas
    vector<Minisat::Lit> relaxation_lits_1;
    for (auto& var : relaxationVars) {
        relaxation_lits_1.push_back(Minisat::mkLit(var));
    }
    encodeRelaxedDeadCells(solver, cellsLits, relaxation_lits_1, k);
    enforceRelaxationBound(solver, relaxationVars, k);

    // bordas são todas mortas
    for (size_t i = 0; i < n; ++i) {
        solver.addClause(~Minisat::mkLit(previous[i][0]));
        solver.addClause(~Minisat::mkLit(previous[i][m - 1]));
    }
    for (size_t j = 0; j < m; ++j) {
        solver.addClause(~Minisat::mkLit(previous[0][j]));
        solver.addClause(~Minisat::mkLit(previous[n - 1][j]));
    }

    for (int i = 1; i < n - 1; ++i) { // Ignorar bordas
        for (int j = 1; j < m - 1; ++j) {
            Minisat::Lit cell = Minisat::mkLit(previous[i][j]);

            // Cria um conjunto com os vizinhos da celula atual
            vector<Minisat::Lit> neighbors;
            for (int di = -1; di <= 1; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                    if (di != 0 || dj != 0) {
                        neighbors.push_back(
                            Minisat::mkLit(previous[i + di][j + dj])
                        );
                    }
                }
            }

            if (board[i][j] == 1) { // vivo
                loneliness(solver, neighbors);
                stagnation(solver, neighbors, cell);
                overcrowding(solver, neighbors);
            } else { // morto
                preservation(solver, neighbors, cell);
                life(solver, neighbors);
            }
        }
    }

    solver.simplify();

    cerr << "Resolvendo... k= " << k << "\n";

    bool isSat = solver.solve();

    if (isSat) {
        vector<vector<int>> result(n, vector<int>(m, 0));

        for (int i = 1; i < n - 1; ++i) { // Ignorar bordas
            for (int j = 1; j < m - 1; ++j) {
                result[i][j] =
                    (solver.modelValue(previous[i][j])
                             == (Minisat::lbool((uint8_t)0))
                         ? 1
                         : 0);
            }
        }

        resultBoards[k] = result;
    }

    return isSat;
}

int main(void) {
    int n, m;
    cin >> n >> m;
    n += 2;
    m += 2; // adiciona as bordas

    std::map<int, vector<vector<int>>> resultBoards;
    addUnsatResult(n, m, resultBoards);

    vector<vector<int>> board(n, vector<int>(m, 0));
    for (int i = 1; i < n - 1; ++i) {
        for (int j = 1; j < m - 1; ++j) {
            cin >> board[i][j];
        }
    }

    int low = 0, high = n * m;
    int best = -1;

    while (low <= high) {
        int k = low + (high - low) / 2;
        if (solve(board, k, resultBoards)) {
            best = k;
            high = k - 1;
        } else {
            low = k + 1;
        }
    }

    if (best != -1) {
        cerr << "valor ideal é " << best << "\n";
    } else {
        cerr << "Nenhuma solução encontrada para qualquer valor de t.\n";
    }

    cout << n - 2 << " " << m - 2 << endl;
    for (int i = 1; i < n - 1; ++i) { // Ignorar bordas
        for (int j = 1; j < m - 1; ++j) {
            cout << resultBoards[best][i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}