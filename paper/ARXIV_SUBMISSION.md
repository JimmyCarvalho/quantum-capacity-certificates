# Submissão ao arXiv: passo a passo

Arquivo para enviar: `physics/paper/arxiv.tar.gz`. Ele contém `main.tex`, `main.bbl`, `refs.bib`,
`thresholds.png` e a pasta `anc/`, com o código, os estados e o script de verificação. O arXiv publica
essa pasta como "ancillary files". A compilação do pacote foi testada: 14 páginas, com as referências.

## 1. Conta
- Crie a conta em https://arxiv.org/user/register com o nome "Jimmy Kelvyn Estevam Carvalho" e a
  afiliação "RCO Tecnologia Ltda.".
- Use um e-mail que você leia sempre, porque é para ele que vão os avisos do arXiv.

## 2. Endosso (endorsement)
- Quem submete pela primeira vez em quant-ph precisa ser "endossado" por um autor da área.
- Na primeira submissão, o arXiv mostra um código de endosso de 6 letras. Mande esse código para um
  pesquisador que já publicou em quant-ph. O endosso não é revisão do artigo; a pessoa só confirma que
  o texto é um trabalho científico sério da área.
- O pedido costuma dar certo quando a pessoa já leu o artigo. Por isso o ideal é que seja o mesmo
  especialista que fizer a revisão. O rascunho do e-mail está no fim deste arquivo.

## 3. Formulário de submissão
- **License:** CC BY 4.0
- **Primary category:** quant-ph (Quantum Physics)
- **Cross-list (opcional):** math-ph
- **Upload:** `arxiv.tar.gz`. Confira o PDF gerado pelo arXiv antes de confirmar.
- **Title** (copiar e colar):

  Emergent Gaussian channels for Pauli noise: certified quantum-capacity thresholds and positive capacity of the additive Gaussian noise channel beyond $1/\sqrt{e}$

- **Authors:** Jimmy Kelvyn Estevam Carvalho
- **Comments:**

  14 pages, 3 figures, 6 tables. Computer-assisted proofs; code, input states and certificate logs are included as ancillary files and at https://github.com/JimmyCarvalho/quantum-capacity-certificates

- **Abstract** (copiar e colar; o arXiv aceita $...$):

  A qubit Pauli channel acting on permutation-symmetric states of $M$ qubits with few excitations around a pole behaves, for large $M$, like a single-mode bosonic Gaussian channel, as Ahmed et al. showed for the depolarizing channel. We give this emergent channel for an arbitrary Pauli channel, show that for the independent X-Z (BB84) channel it is the additive Gaussian noise (AGN) channel, and show that it is antidegradable exactly when the transverse shrinking factor is at most $1/\sqrt{2}$, the phase-covariant cloning limit. With an exact Dicke-basis evolution and interval-arithmetic certifiers we certify positive coherent information of explicit states: BB84 at $q=0.12675$ (240 qubits) and the 2-Pauli channel at $p=0.1255$ (240 qubits), above the previous $0.118371$ and $0.118067$, and the depolarizing channel at total error $0.1977936$ (240 qubits); in the Gaussian limit we certify $0.1997875$, the largest depolarizing bound with a certified final sign that we know of (the larger $0.2029161$ rests on a floating-point value). In the Gaussian limit we certify $Q>0$ for the AGN channel for all noise variances $\sigma^2\le 0.4121$ ($\sigma\le0.6419$), beyond $\sigma^2=1/e$ where Gaussian inputs stop; this refutes the 2001 speculation of Harrington and Preskill that the Gaussian-input coherent information is the capacity of this channel. More than half of this gain already comes from a three-parameter squeezed cat. For the thermal attenuator we certify $Q>0$ inside the open window between the Gaussian-input threshold and antidegradability at every temperature tried, closing 24-50% of it (e.g. $\eta\ge0.77875$ with one thermal photon, below the recently certified $0.7841$), and we give a thermal amplifier with $Q>0$ where all Gaussian inputs fail. All results are computer-assisted; the certificates were re-computed with independent code.

## 4. Depois de enviar
- O artigo aparece no arXiv no dia útil seguinte, por volta das 21h (horário de Brasília). A
  moderação pode reter ou reclassificar o artigo.
- Quando sair o número (arXiv:2610.xxxxx), troque o link do tweet pelo link do arXiv
  (https://arxiv.org/abs/2610.xxxxx), que tem mais credibilidade que o do GitHub.
- Para atualizar o artigo depois (por exemplo, depois da revisão), use "Replace", que gera a versão v2.

## Rascunho de e-mail para pedir revisão e endosso (em inglês)

Subject: Request for a quick look and arXiv endorsement (quant-ph): Q > 0 for the AGN channel beyond 1/sqrt(e)

Dear Prof. [Name],

I am an independent researcher and would like to submit to quant-ph a computer-assisted paper whose main
results are (i) positive quantum capacity of the additive Gaussian noise channel for all
sigma^2 <= 0.4121, beyond the Holevo-Werner threshold 1/e, (ii) Q > 0 for the thermal attenuator with one
thermal photon down to eta = 0.77875, and (iii) new lower bounds on the quantum-capacity thresholds of the
BB84 and 2-Pauli channels. Every claim is an interval-arithmetic certificate (python-flint/Arb) for an
explicit input state; code, states and logs are public:
https://github.com/JimmyCarvalho/quantum-capacity-certificates
(paper: paper/main.pdf). The computations were done with the help of an AI system, and the paper says so.

Since this is my first arXiv submission, I need an endorsement for quant-ph. If you find the work sound
after a look, would you be willing to endorse it (code: [XXXXXX])? Any comments or problems you spot would
be very welcome.

Best regards,
Jimmy Kelvyn Estevam Carvalho
RCO Tecnologia Ltda.
