#!/bin/tcsh

set n = 0

foreach file (`find . -type f`)
    git add "$file"
    @ n++

    if ($n % 1000 == 0) then
        echo "Committing batch of $n files..."
        git commit -m "Pushing batch of $n files"
        git push origin main
    endif
end

if ($n % 1000 != 0) then
    echo "Committing final batch of $n files..."
    git commit -m "Pushing final batch of $n files"
    git push origin main
endif

