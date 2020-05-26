header = """\\captionsetup[figure]{font=footnotesize,labelfont=footnotesize}
\\begin{figure}
\\begin{tabular}{cccccc}
"""

body ="""
\\resizebox{0.125\\textwidth}{!}{
\\begin{tikzpicture}
	\\node[circle,draw=black,fill=black!20] (1) at (2,0) {\\textbf{1}};
	\\node[circle,draw=black,fill=black!20] (2) at (1,1.732) {\\textbf{2}};
	\\node[circle,draw=black,fill=black!20] (3) at (-1,1.732) {\\textbf{3}};
	\\node[circle,draw=black,fill=black!20] (4) at (-2,0) {\\textbf{4}};
	\\node[circle,draw=black,fill=black!20] (5) at (-1,-1.732) {\\textbf{5}};
	\\node[circle,draw=black,fill=black!20] (6) at (1,-1.732) {\\textbf{6}};
	\\foreach \\x/\\y in {1/4,1/3}
		\\draw[black,=>latex',->, thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {2/3,2/4}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {3/1,3/2,3/4,3/5}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {4/1,4/2,4/3,4/6}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {5/1,5/2,5/4,5/6}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {6/1,6/2,6/3,6/5}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
\\end{tikzpicture}%
} 
& 
 \\includegraphics[width=0.2\\linewidth]{images/case3_1.png}
&
\\resizebox{.125\\textwidth}{!}{
\\begin{tikzpicture}
	\\node[circle,draw=black,fill=black!20] (1) at (2,0) {\\textbf{1}};
	\\node[circle,draw=black,fill=black!20] (2) at (1,1.732) {\\textbf{2}};
	\\node[circle,draw=black,fill=black!20] (3) at (-1,1.732) {\\textbf{3}};
	\\node[circle,draw=black,fill=black!20] (4) at (-2,0) {\\textbf{4}};
	\\node[circle,draw=black,fill=black!20] (5) at (-1,-1.732) {\\textbf{5}};
	\\node[circle,draw=black,fill=black!20] (6) at (1,-1.732) {\\textbf{6}};
	\\foreach \\x/\\y in {1/4,1/3}
		\\draw[black,=>latex',->, thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {2/3,2/4}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {3/1,3/2,3/4,3/5}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {4/1,4/2,4/3,4/6}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {5/1,5/2,5/4,5/6}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {6/1,6/2,6/3,6/5}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
\\end{tikzpicture}%
} 
& 
 \\includegraphics[width=0.2\\linewidth]{images/case3_1.png}
&
\\resizebox{.125\\textwidth}{!}{
\\begin{tikzpicture}
	\\node[circle,draw=black,fill=black!20] (1) at (2,0) {\\textbf{1}};
	\\node[circle,draw=black,fill=black!20] (2) at (1,1.732) {\\textbf{2}};
	\\node[circle,draw=black,fill=black!20] (3) at (-1,1.732) {\\textbf{3}};
	\\node[circle,draw=black,fill=black!20] (4) at (-2,0) {\\textbf{4}};
	\\node[circle,draw=black,fill=black!20] (5) at (-1,-1.732) {\\textbf{5}};
	\\node[circle,draw=black,fill=black!20] (6) at (1,-1.732) {\\textbf{6}};
	\\foreach \\x/\\y in {1/4,1/3}
		\\draw[black,=>latex',->, thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {2/3,2/4}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {3/1,3/2,3/4,3/5}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {4/1,4/2,4/3,4/6}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {5/1,5/2,5/4,5/6}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
	\\foreach \\x/\\y in {6/1,6/2,6/3,6/5}
		\\draw[black,=>latex',->,thick] (\\x) -- (\\y);
\\end{tikzpicture}%
} 
& 
 \\includegraphics[width=0.2\\linewidth]{images/case3_1.png}\\\\ 
 
\\scriptsize{$G_{P.1}$}  & \\scriptsize{$W_r(L(G_{P.2}))$}&\\scriptsize{$G_{P.1}$}  
& \\scriptsize{$W_r(L(G_{P.2}))$} &\\scriptsize{$G_{P.1}$}  & \\scriptsize{$W_r(L(G_{P.2}))$} \\\\[6pt]
"""



footer =  """\\end{tabular}
\\end{figure}"""


#can fit 8 per page

outfile = open('latex_outfile.txt', 'w+')
for j in range(2):
    outfile.write(header)
    for i in range(8):
        outfile.write(body)
    outfile.write(footer)


