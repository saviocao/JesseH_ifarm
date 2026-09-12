#!/bin/tcsh

foreach gitdir (`find . -type d -name .git -not -path './.git'`)
    echo "Removing: $gitdir"
    rm -rf "$gitdir"
end

echo "Done."

















