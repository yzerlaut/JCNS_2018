
def find_boundaries(filename='paper.tex', folder='tex/'):
    f = open(folder+filename)
    l = f.readline()
    i, i_limit = 0, 10000

    while (l!='\\section{Abstract}\n') and (l!='\\section{Key points summary}\n') and (i<i_limit):
        l = f.readline()
        i+=1
    MANUSCRIPT_START = i
    while (l!='\\section{Introduction}\n') and (i<i_limit):
        l = f.readline()
        i+=1
    INTRODUCTION_START = i+1
    while (l!='\\section{References}\n') and (l!='\\begin{thebibliography}{}\n') and (i<i_limit):
        l = f.readline()
        i+=1
    MANUSCRIPT_END = i
    return MANUSCRIPT_START, INTRODUCTION_START, MANUSCRIPT_END

def reformat_line(line):
    if len(line.split('\\label{sec-'))>1:
        return ''
    else:
        return line

def produce_tex_file(filename='paper', folder='tex/', replace=True,\
                     with_biblio=False, full_file=False):
    MANUSCRIPT_START, INTRODUCTION_START, MANUSCRIPT_END = \
      find_boundaries(filename=filename, folder=folder)
    f = open(folder+filename)
    # empty read
    for i in range(MANUSCRIPT_START):
        f.readline()
    # manuscript read
    core_manuscript = ''

    if full_file:
        core_manuscript += '\\documentclass[a4paper, colorlinks]{article} \n'
        core_manuscript += '\\usepackage{hyperref, lineno, graphicx} \n'
        core_manuscript += '\\usepackage[utf8]{inputenc} \n'
        core_manuscript += '\\renewcommand{\\includegraphics}[2][]{\\fbox{#2}} % remove graphics \n'
        core_manuscript += '\\begin{document} \n'

    for i in range(INTRODUCTION_START-MANUSCRIPT_START):
        core_manuscript += reformat_line(f.readline())
    core_manuscript += '\\linenumbers \n'
    for i in range(MANUSCRIPT_END-INTRODUCTION_START):
        core_manuscript += reformat_line(f.readline())
    core_manuscript += '\\nolinenumbers \n'

    if with_biblio:
        core_manuscript += open(folder+filename+'.bbl').read()
    else:
        core_manuscript += '\\bibliographystyle{apalike}\n'
        core_manuscript += '\\bibliography{biblio}\n'
    
    if full_file:
        core_manuscript += '\\end{document} \n'
    
    # then some replacements
    if replace:
        core_manuscript = core_manuscript.replace('./figures/', '../figures/')
        # core_manuscript = core_manuscript.replace('.png', '.eps')
        core_manuscript = core_manuscript.replace('\\section', '\\section*')
        core_manuscript = core_manuscript.replace('\\subsection', '\\subsection*')
        core_manuscript = core_manuscript.replace('\\citetext', '\\cite')
        core_manuscript = core_manuscript.replace('\\textcolor{red}', '\\textbf')

    new_paper = open(folder+'simple_paper.tex', 'w')
    new_paper.write(core_manuscript)
    # new_paper.write("\n \nolinenumbers")
    new_paper.close()

def run_compilation(filename='paper', folder='tex/'):
    
    import os
    
    os.system('pdflatex -shell-escape -interaction=nonstopmode '+folder+filename+'.tex')
    os.system('bibtex -terse '+folder+filename+'.aux')
    os.system('pdflatex -shell-escape -interaction=nonstopmode '+folder+filename+'.tex')
    os.system('pdflatex -shell-escape -interaction=nonstopmode '+folder+filename+'.tex')
    
    return None

if __name__=='__main__':

    produce_tex_file(filename='paper.tex', folder='./',\
                     with_biblio=False, full_file=True, replace=True)
    run_compilation(filename='simple_paper', folder='./')
