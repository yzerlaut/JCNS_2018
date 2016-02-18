
def find_boundaries(filename='paper.tex', folder='tex/'):
    f = open(folder+filename)
    l = f.readline()
    i, i_limit = 0, 10000
    while (i<i_limit):
        l = f.readline()
        if l=='\\section{Abstract}\n':
            MANUSCRIPT_START = i+1
        if l=='\\section{Introduction}\n':
            INTRODUCTION_START = i+1
            print i
        if l=='\\section{References}\n' or l=='\\begin{thebibliography}{}\n':
            MANUSCRIPT_END = i
            i=i_limit
        i+=1
    return MANUSCRIPT_START, INTRODUCTION_START, MANUSCRIPT_END

def reformat_line(line):
    if len(line.split('\\label{sec-'))>1:
        return ''
    else:
        return line

def produce_tex_file(filename='paper.tex', folder='tex/', replace=True,\
                     with_biblio=False):
    MANUSCRIPT_START, INTRODUCTION_START, MANUSCRIPT_END = \
      find_boundaries(filename=filename, folder=folder)
    f = open(folder+filename)
    # empty read
    for i in range(MANUSCRIPT_START):
        f.readline()
    # manuscript read
    core_manuscript = ''
    for i in range(INTRODUCTION_START-MANUSCRIPT_START):
        core_manuscript += reformat_line(f.readline())
    core_manuscript += '\\linenumbers \n'
    for i in range(MANUSCRIPT_END-INTRODUCTION_START):
        core_manuscript += reformat_line(f.readline())
    core_manuscript += '\\nolinenumbers \n'

    if with_biblio:
        core_manuscript += open(folder+'paper.bbl').read()
    else:
        core_manuscript += '\\bibliography{tex/biblio}\n'
    
    core_manuscript += '\\end{document} \n'
    
    # then some replacements
    if replace:
        core_manuscript = core_manuscript.replace('./figures/', '../figures/')
        core_manuscript = core_manuscript.replace('.png', '.eps')
        core_manuscript = core_manuscript.replace('\\section', '\\section*')
        core_manuscript = core_manuscript.replace('\\subsection', '\\subsection*')
        core_manuscript = core_manuscript.replace('\\citetext', '\\cite')
        core_manuscript = core_manuscript.replace('\\textcolor{red}', '\\textbf')

    new_paper = open(folder+'simple_paper.tex', 'w')
    new_paper.write(core_manuscript)
    # new_paper.write("\n \nolinenumbers")
    new_paper.close()

def run_compilation(filename='paper.tex', folder='tex/'):
    
    import os
    
    produce_tex_file(filename=filename)
    os.system('latex -shell-escape -interaction=nonstopmode -output-directory=tex '+folder+'new_paper.tex')
    os.system('bibtex -terse '+folder+'new_paper.aux')
    produce_tex_file(filename=filename, with_biblio=True, replace=False)
    os.system('latex -shell-escape -interaction=nonstopmode -output-directory=tex '+folder+'new_paper.tex')
    os.system('bibtex -terse '+folder+'new_paper.aux')
    os.system('latex -shell-escape -interaction=nonstopmode -output-directory=tex '+folder+'new_paper.tex')
    os.system('latex -shell-escape -interaction=nonstopmode -output-directory=tex '+folder+'new_paper.tex')
    os.system('dvipdf '+folder+'new_paper.dvi')
    os.system('mv new_paper.pdf pdf_output/paper.pdf')
    
    return None
