# Time-Reversal in Conway’s Life as SAT, while minimizing living cells

Este projeto explora o jogo da vida de Conway, um autômato celular famoso onde células em uma grade vivem, morrem ou nascem com base em regras simples de vizinhança. 
A ideia aqui é inverter o tempo: dado um estado, tentamos descobrir um estado anterior possível que evolui para ele — usando SAT (satisfatibilidade booleana). 
Além disso, buscamos uma solução com o menor número possível de células vivas no início.

Essa tarefa exige muito esforço computacional. Dependendo do estado dado, meu projeto consegue solucionar entradas de tamanho até 20x20 em poucos minutos.
