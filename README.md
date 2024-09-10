# LSound

Código acadêmico para a simulação de problemas de acústica linear utilizando o método dos elementos finitos. 


## Informações Gerais

O pré e o pós processamento são realizados utilizando o gmsh (https://gmsh.info/).

Os seguintes elementos são disponíveis: 

+ bilinear isoparamétrico de 4 nós

+ triangular de 3 nós

+ trilinear hexaédrico de 8 nós 

+ tetraedro linear de 4 nós. 

Materiais e condição de contorno são informadas por meio de Physical Groups do gmsh. O programa atualmente reconhece as seguintes informações:

> Material, $id$, $\rho$, $c$, $Z$

Informa as propriedades do material. $id$ é um inteiro utilizado para identificar o material, $\rho$ a densidade, $c$ a velocidade do som e $Z$ 

> Open

Informa nós com pressão nula (condição aberto)

> Vn, A, freq, fase 

Informa fronteira com velocidade normal imposta (pistão). A velocidade é assumida na forma $v_n(t) = A\cos(freq*t + fase)$ onde a frequência é em Hz e a fase é informada em graus.

> Yn,valor

Informa o amortecimento.

> Monitor 

Informa nós para serem exportados na análise harmônica.

## Tipos de análise

A rotina principal do programa é 

> Analise(meshfile,método)

onde $meshfile$ é um arquivo de entrada $.msh$ do gmsh e método uma das seguintes opções:

+ :Modal

+ :Harmonic

+ :Newmark

+ :Bathe

Cada método tem um conjunto de opções associado:

> $f$, $X$ = Analise(meshfile,:Modal,nev)

onde $nev$ é o número de autovalores a serem calculados, $f$ um vetor com as frequências naturais (em Hz) e $X$ uma matriz com os autovetores. Essas informações são exportadas para um arquivo chamado $Modal.pos$


> monitor, $X$ = Analise(meshfile,:Harmonic,freqs)

onde $freqs$ é uma lista com frequências a serem analizadas (em Hz). $monitor$ é um vetor com os nós que estão sendo monitorados e $X$ é uma matriz complexa com os nós como linhas e frequências como colunas. 

> $T$,$X$ = Analise(meshfile,:Newmark,Tf,$\Delta t$, $U0$, $V0$, output)

onde $Tf$ é o tempo final, $\Delta t$ a discretização no tempo, $U0$ e $V0$ as condições iniciais e $output$ um flag lógico para gravar ou não gravar as informações no gmsh (o arquivo pode ser bem grande...). $T$ é um vetor com os tempos discretos que foram analisados e $X$ uma matriz com os nós (linhas) em cada tempo (colunas). 

> $T$,$X$ = Analise(meshfile,:Bathe,Tf,$\Delta t$, $U0$, $V0$, output)

onde $Tf$ é o tempo final, $\Delta t$ a discretização no tempo, $U0$ e $V0$ as condições iniciais e $output$ um flag lógico para gravar ou não gravar as informações no gmsh (o arquivo pode ser bem grande...). $T$ é um vetor com os tempos discretos que foram analisados e $X$ uma matriz com os nós (linhas) em cada tempo (colunas). 