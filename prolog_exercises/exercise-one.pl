% some facts:

student(mary).
student(christine).
student(john).
student(peter).

teacher(rob).
teacher(lidia).

%---- All teachers like first-order logic (l)

likes(X, logic) :- teacher(X).

%---- Some student likes first-order logic (l)

% we can't express the abstract concept "some" in prolog, we need
% concrete examples:

likes(john, logic).
likes(mary, logic).

% ask prolog "which students like logic?"
% this is the same as asking:
% "does there exist at least one X such that student(X) and likes(X,logic)?"

% ?- likes(X, logic), student(X).

% who likes logic?

% ?- likes(X, logic). % can try alternative solutions until fail (false)

% Does John like logic?

% ?- likes(john, logic). % YES

%---- no student is happy = all students are unhappy

unhappy(X) :- student(X).
happy(X) :- not(unhappy(X)).

% who is not happy?
