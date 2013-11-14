#!/usr/bin/env python

import sys

current_word = None
count_curr = 0
word = None

for line in sys.stdin:
    #Remove leading and trailing whitespace
    line = line.strip()

    word, count = line.split('\t',1)
    try:
        count = int(count)
    except ValueError:
        continue
    if current_word == word:
        count_curr += count
    else:
        if current_word:
            print '%s,%i' % (current_word, count_curr)
        count_curr = count
        current_word = word

if current_word == word:
    print '%s,%i' % (current_word, count_curr)

exit()