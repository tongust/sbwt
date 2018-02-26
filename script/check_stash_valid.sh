git status | grep "mod" |grep "[A-Za-z_]*\.[a-z_]*" -o | xargs cat | md5sum
