#!/usr/bin/env bash

# Text colors
if [[ $TERM != *"xterm"* ]]; then
    export TERM=dummy
fi

# colors for the foreground texts
fgred=$(tput setaf 1)    # Red
fggreen=$(tput setaf 2)    # Green
fgyellow=$(tput setaf 3)    # Yellow
fgblue=$(tput setaf 4)    # Blue
fgpurple=$(tput setaf 5)    # Purple
fgcyan=$(tput setaf 6)    # Cyan
fgwhite=$(tput setaf 7)    # White
fgrst=$(tput sgr0)       # Text reset

# colors for the background text
bgred=$(tput setab 1)    # Red
bggreen=$(tput setab 2)    # Green
bgyellow=$(tput setab 3)    # Yellow
bgblue=$(tput setab 4)    # Blue
bgpurple=$(tput setab 5)    # Purple
bgcyan=$(tput setab 6)    # Cyan
bgwhite=$(tput setab 7)    # White
bgrst=$(tput sgr0)       # reset
