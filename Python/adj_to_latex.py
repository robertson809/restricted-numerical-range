import numpy as np
from polyGraphs import plt_nr
n = 6
mat_adj_file= open('true_poly_6_adj.txt')
mat_lap_file = open('true_poly_6_lap.txt')
# outfile = open('true_poly_6_lap.txt', 'w+')
adj_lineList = mat_adj_file.readlines()
lap_lineList = mat_lap_file.readlines()
adj = np.zeros((n,n),dtype=float)
lap = np.zeros((n,n),dtype=float)
i = 0
count = 0
for_each_str_list = []
# for adj_row in adj_lineList:
for adj_row, lap_row in zip(adj_lineList, lap_lineList):
    for_each_str = ''
    # lap_row = adj_row
    adj_row = adj_row.split(" ")
    lap_row = lap_row.split(" ")
    if len(adj_row) < n:
        count += 1
        print('Building graph {} of 111'.format(count))
        # write the latex graphing line
        for row_i in range(n):
            for column_j in range(n):
                if int(adj[row_i][column_j]) == 1:
                    for_each_str += '{}/{},'.format(row_i+1, column_j+1)
        for_each_str_list.append(for_each_str[:-1]) # -1 is to remove the last comma

        # create the figure
        # plt_nr(lap, save_num = count)

        # used for converting adj to lap
        # diag = np.array([np.sum(adj[i, :]) for i in range(n)])
        # lap = np.diag(diag) - adj
        # for k in range(n):
        #     for j in range(n - 1):
        #         outfile.write("%d " % lap[k, j])
        #     outfile.write("%d\n" % lap[k, n - 1])
        # outfile.write("\n")


        i = 0
    else:
        # store row
        adj[i, :] = [eval(adj_row[j]) for j in range(n)]
        lap[i, :] = [eval(lap_row[j]) for j in range(n)]
        # update i
        i = i + 1
# outfile.close()
header = """\\captionsetup[figure]{font=footnotesize,labelfont=footnotesize}
\\begin{figure}
\\begin{tabular}{cccccc}
"""

body ="""
\\resizebox{0.1\\textwidth}{!}{
\\begin{tikzpicture}
	\\node[circle,draw=black,fill=black!20] (1) at (2,0) {\\textbf{1}};
	\\node[circle,draw=black,fill=black!20] (2) at (1,1.732) {\\textbf{2}};
	\\node[circle,draw=black,fill=black!20] (3) at (-1,1.732) {\\textbf{3}};
	\\node[circle,draw=black,fill=black!20] (4) at (-2,0) {\\textbf{4}};
	\\node[circle,draw=black,fill=black!20] (5) at (-1,-1.732) {\\textbf{5}};
	\\node[circle,draw=black,fill=black!20] (6) at (1,-1.732) {\\textbf{6}};
	\\foreach \\x/\\y in {%s}
		\\draw[-{Latex[length=3mm]}] (\\x) -- (\\y);
\\end{tikzpicture}
} 
& 
 \\includegraphics[width=0.1\\linewidth]{appendix_polygraphs/polyGraph%s.png}
&
\\resizebox{.1\\textwidth}{!}{
\\begin{tikzpicture}
	\\node[circle,draw=black,fill=black!20] (1) at (2,0) {\\textbf{1}};
	\\node[circle,draw=black,fill=black!20] (2) at (1,1.732) {\\textbf{2}};
	\\node[circle,draw=black,fill=black!20] (3) at (-1,1.732) {\\textbf{3}};
	\\node[circle,draw=black,fill=black!20] (4) at (-2,0) {\\textbf{4}};
	\\node[circle,draw=black,fill=black!20] (5) at (-1,-1.732) {\\textbf{5}};
	\\node[circle,draw=black,fill=black!20] (6) at (1,-1.732) {\\textbf{6}};
	\\foreach \\x/\\y in {%s}
		\\draw[-{Latex[length=3mm]}] (\\x) -- (\\y);
\\end{tikzpicture}
} 
& 
 \\includegraphics[width=0.1\\linewidth]{appendix_polygraphs/polyGraph%s.png}
&
\\resizebox{.1\\textwidth}{!}{
\\begin{tikzpicture}
	\\node[circle,draw=black,fill=black!20] (1) at (2,0) {\\textbf{1}};
	\\node[circle,draw=black,fill=black!20] (2) at (1,1.732) {\\textbf{2}};
	\\node[circle,draw=black,fill=black!20] (3) at (-1,1.732) {\\textbf{3}};
	\\node[circle,draw=black,fill=black!20] (4) at (-2,0) {\\textbf{4}};
	\\node[circle,draw=black,fill=black!20] (5) at (-1,-1.732) {\\textbf{5}};
	\\node[circle,draw=black,fill=black!20] (6) at (1,-1.732) {\\textbf{6}};
	\\foreach \\x/\\y in {%s}
		\\draw[-{Latex[length=3mm]}] (\\x) -- (\\y);
\\end{tikzpicture}
} 
& 
 \\includegraphics[width=0.1\\linewidth]{appendix_polygraphs/polyGraph%s.png}\\\\ 
 
\\scriptsize{$G_{%s}$}  & \\scriptsize{$W_r(L(G_{%s}))$}&\\scriptsize{$G_{%s}$}  
& \\scriptsize{$W_r(L(G_{%s}))$} &\\scriptsize{$G_{%s}$}  & \\scriptsize{$W_r(L(G_{%s}))$} \\\\[6pt]
"""



footer =  """\\end{tabular}
\\end{figure}"""


#can fit 8 per page

outfile = open('latex_outfile.txt', 'w+')
graph_count = 0
for j in range(4):
    outfile.write(header)
    for i in range(8):
        outfile.write(body % (for_each_str_list[graph_count], graph_count +1, for_each_str_list[graph_count+1],
                              graph_count+2, for_each_str_list[graph_count+2], graph_count+3,
                              graph_count + 1, graph_count + 1, graph_count + 2,
                              graph_count + 2, graph_count + 3, graph_count + 3
                              ))
        graph_count += 3
    outfile.write(footer)
#remainder
outfile.write(header)
for i in range(5):
    outfile.write(body % (for_each_str_list[graph_count], graph_count + 1, for_each_str_list[graph_count+1],
                              graph_count+2, for_each_str_list[graph_count+2], graph_count + 3,
                          graph_count + 1, graph_count+2, graph_count + 3,
                          graph_count + 1, graph_count+2, graph_count + 3))
    graph_count += 3
print("final graph count is", graph_count)
outfile.write(footer)


