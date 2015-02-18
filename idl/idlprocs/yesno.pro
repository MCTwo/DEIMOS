function yesno, prompt=prompt
;+
; NAME:
;      yesno
;
;
; PURPOSE: 
;      This function reads the response to a yes or no question and
;      returns a boolean T/F
;
; CALLING SEQUENCE: 
;      result=yesno([prompt=])
;
;
; KEYWORD PARAMETERS: 
;      prompt:  The question that is to be asked, if the default string,
;               'Yes or no?', is not desired.
;
;
; OUTPUTS:
;      result:
;      If response to prompt is a variation of the word yes (including
;         y,  Yes, YEs, etc), returns 1
;      If response to prompt is a variation of the word no (including
;         n, NO, No, etc.), returns 0
;      If response to prompt is neither of these, returns -1
;
;
; EXAMPLE:
;      IDL> print, yesno()
;      Yes or no?
;      : y
;            1
; EXAMPLE:
;      IDL> print, yesno(prompt='Are you cool?')
;      Are you cool?
;      : NO
;            0
;
;
; MODIFICATION HISTORY: 
;      Written by: Cassandra VanOutryve, July 2003
;
;-



if keyword_set(prompt) then print, prompt else print, 'Yes or no?'

a=' '
read, a; read the response

b=strlowcase(a); put the response into all lower case letters
c=strmid(b, 0, 1); take just the first character of the string

case c of
    'y': return, 1B
    'n': return, 0B
    else: return, -1
endcase

end

