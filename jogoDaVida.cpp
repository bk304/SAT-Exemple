#include <bits/stdc++.h>
using namespace std;

/*
As regras que definem a evolução de um estado para outro do jogo são:
  - toda célula morta com exatamente três vizinhos vivos torna-se viva;
  - toda célula viva com menos de dois vizinhos vivos morre;
  - toda célula viva com mais de três vizinhos vivos morre;
  - toda célula viva com dois ou três vizinhos vivos permanece viva.
*/

#define MORTO 0
#define VIVO 1

int main(void){
    int n,m;
    cin >> n >> m;
    n+=2; m+=2; // adiciona as bordas

    vector<vector<int>> board_t0(n, vector<int>(m, 0));
    vector<vector<int>> board_t1(n, vector<int>(m, 0));
    for(int i=1; i<n-1; ++i){
        for(int j=1; j<m-1; ++j){
            cin >> board_t0[i][j];
        }
    }

    for (int i = 1; i < n - 1; ++i) { // Ignorar bordas
        for (int j = 1; j < m - 1; ++j) {

            int neighbors = 0;
            for (int di = -1; di <= 1; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                    if (di != 0 || dj != 0) {
                        if(board_t0[i+di][j+dj] == VIVO)
                            neighbors += 1;
                    }
                }
            }

            if(board_t0[i][j] == MORTO && neighbors == 3)
                board_t1[i][j] = VIVO;
            
            if(board_t0[i][j] == VIVO && neighbors < 2)
                board_t1[i][j] = MORTO;

            if(board_t0[i][j] == VIVO && neighbors > 3)
                board_t1[i][j] = MORTO;

            if(board_t0[i][j] == VIVO && (neighbors == 2 || neighbors == 3))
                board_t1[i][j] = VIVO;
        }
    }

    cout << n -2 << " " << m -2 << endl;
    for (int i = 1; i < n - 1; ++i) { // Ignorar bordas
        for (int j = 1; j < m - 1; ++j) {
            cout << (board_t1[i][j] == VIVO ? 1 : 0) << " ";
        }
        cout << endl;
    }

    return 0;
}