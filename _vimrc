cd d:\lrcode\cpp
set guifont=Consolas:h12:cANSI

set guioptions-=m
set guioptions-=T

source $VIMRUNTIME/mswin.vim
behave mswin

set nu
set softtabstop=4
set tabstop=4
set shiftwidth=4
set nobackup
set mouse=a
set hlsearch
set ruler

map <F2> : call Save() <CR>
map <C-x> : !mvim -ro % <CR>
map <F5> : !gdb --quiet %< <CR>
map <F9> : call Compile() <CR>
map <F8> : !bash check.sh <CR>
map <F7> : !%<.exe <.in <CR>
map <F10> : call Run() <CR>
map <F4> : !echo % <CR>

syn keyword Type dint
syn keyword Special PROC
syn keyword Special lld
filetype plugin indent on
syntax enable

colors koehler

"set makeprg=g++\ %\ -o\ %<\ -g\ -Wall\ -Wextra\ -DLAEKOV_LOCAL

func Save()
	if &filetype=='javascript'
		set expandtab
	endif
	exec "w"
endfunc

func Compile()
	exec "w"
	if &filetype=='cpp'
		exec "!g++ % -o %< -g -Wall -Wextra -DLAEKOV_LOCAL && size %<.exe"
	elseif &filetype=='c'
		exec "!gcc % -o %< -g"
	elseif &filetype=='pascal'
		exec "!fpc % -g"
	elseif &filetype=='sh'
		exec "!chmod +x %"
	elseif &filetype=='java'
		exec "!javac %"
	elseif &filetype=='tex'
		exec "!xelatex %"
	endif
endfunc

func Run()
	if &filetype=='python'
		exec "!python %"
	elseif &filetype=='sh'
		exec "!bash %"
	elseif &filetype=='java'
		exec "!java Main"
	elseif &filetype=='html'
		exec "!open ./%"
	elseif &filetype=='tex'
		exec "!open %<.pdf"
	elseif &filetype=='javascript'
		exec "!node %"
	else
		exec "!%<.exe"
	endif
endfunc

