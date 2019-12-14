# ~/.bashrc: executed by bash(1) for non-login shells.
# see /usr/share/doc/bash/examples/startup-files (in the package bash-doc)
# for examples

# If not running interactively, don't do anything
[ -z "$PS1" ] && return

# don't put duplicate lines in the history. See bash(1) for more options
# ... or force ignoredups and ignorespace
HISTCONTROL=ignoredups:ignorespace

# append to the history file, don't overwrite it
shopt -s histappend

# for setting history len

# set variable identifying the chroot you work in (used in the prompt below)
if [ -z "$debian_chroot" ] && [ -r /etc/debian_chroot ]; then
    debian_chroot=$(cat /etc/debian_chroot)
fi

# set a fancy prompt (non-color, unless we know we "want" color)
case "$TERM" in
    xterm-color) color_prompt=yes;;
esac

# uncomment for a colored prompt, if the terminal has the capability; turned
# off by default to not distract the user: the focus in a terminal window
# should be on the output of commands, not on the prompt
#force_color_prompt=yes

if [ -n "$force_color_prompt" ]; then
    if [ -x /usr/bin/tput ] && tput setaf 1 >&/dev/null; then
	# We have color support; assume it's compliant with Ecma-48
	# (ISO/IEC-6429). (Lack of such support is extremely rare, and such
	# a case would tend to support setf rather than setaf.)
	color_prompt=yes
    else
	color_prompt=
    fi
fi

if [ "$color_prompt" = yes ]; then
    PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
else
    PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
fi
unset color_prompt force_color_prompt

# If this is an xterm set the title to user@host:dir
case "$TERM" in
xterm*|rxvt*)
    PS1="\[\e]0;${debian_chroot:+($debian_chroot)}\u@\h: \w\a\]$PS1"
    ;;
*)
    ;;
esac

# enable color support of ls and also add handy aliases
if [ -x /usr/bin/dircolors ]; then
    test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"
    alias ls='ls --color=auto'
    #alias dir='dir --color=auto'
    #alias vdir='vdir --color=auto'

    alias grep='grep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias egrep='egrep --color=auto'
fi

# some more ls aliases
alias ll='ls -alF'
alias la='ls -A'
alias l='ls -CF'

# Alias definitions.
# You may want to put all your additions into a separate file like
# ~/.bash_aliases, instead of adding them here directly.
# See /usr/share/doc/bash-doc/examples in the bash-doc package.

if [ -f ~/.bash_aliases ]; then
    . ~/.bash_aliases
fi

# enable programmable completion features (you don't need to enable
# this, if it's already enabled in /etc/bash.bashrc and /etc/profile
# sources /etc/bash.bashrc).
#if [ -f /etc/bash_completion ] && ! shopt -oq posix; then
#    . /etc/bash_completion
#fi

#PS1="\[\033[01;34;1m\][\A]\e[36;1m\u@\e[01;32;1m\h[\W$(__git_ps1 " (%s)")]\$ \[\033[0m\] "
export PS1='\[\033[01;32;1m\][\A]\u@\h[\[\033[01;36;1m\]\W$(__git_ps1 " (%s)")\[\033[01;32;1m\]]\$ \[\033[0m\]'

export sub='/mnt/c/Subsys'

#export PS1="\e[36;1m\u@\e[32;1m\H>\e[0m"
#export PS1='[\u@\h \W$(__git_ps1 " (%s)")]\$ '

source /mnt/c/Subsys/software/root/bin/thisroot.sh
source $sub/software/geant4_install/bin/geant4.sh
function mkcd()
{  
mkdir -p $1
cd $1
}

alias mkg4='mkcd build && cmake ../ && make -j8'
#alias cdg4='cd $G4WORKDIR'
alias uenv='source ~/.bashrc'
alias senv='vim ~/.bashrc'
alias cls='clear && ll'

#alias cq='condor_q'
#alias chis='condor_history -match 10'
#alias cs='condor_submit'

alias root='root -l'
alias rt='root -b'
alias rtc='root-config'

export LZW='liziwei@210.45.78.123:/home/liziwei/R710'
export CRYHOME=$sub/software/cry_v1.7
export CADmeshPATH=$sub/software/CADMesh_install
export code=$sub/work/g4code
export gwk=$sub/work/FTOF_lzw
export wave=$sub/OSC_program
export expdata=/mnt/c/Experiment/Data/labtest
export oned=/mnt/c/Users/liziwei/OneDrive

alias cdgk='cd $gwk'
alias grun='./CRTest mac/quartz.gdml mac/test.mac'
alias gr='./CRTest mac/quartz.gdml'
alias gs='git status'
alias gsh='git show'

export DISPLAY=:0.0


