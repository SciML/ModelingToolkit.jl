%% juliacon.cls - version 1.0

%% Inspired by the template from the International Journal of Computer Applications (IJCA)

\usepackage[scaled=0.92]{helvet}
\def\fileversion{v1.0}
\def\filedate{2019 04 07}
%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{juliacon}
\RequirePackage{latexsym}
\RequirePackage{url}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}

\newif\ifmanuscript
\@twosidetrue\@mparswitchtrue
%
\newdimen\trimheight
\newdimen\trimwidth
\newdimen\typeheight
\newdimen\typewidth
\newdimen\normaltextheight
\newdimen\blindfoliodrop
\newbox\tempbox
%%

\input{journal_dat}

%
\frenchspacing % oh lala bravo quelle belle idée
\DeclareOption{manuscript}{\manuscripttrue}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11.69in}%
    \setlength\paperwidth  {8.27in}%
    \def\special@paper{8.5in,11in}
	\special{papersize=8.5in,11in}}

\DeclareOption{openbib}{%
  \AtEndOfPackage{%
   \renewcommand\@openbib@code{%
      \advance\leftmargin\bibindent
      \itemindent -\bibindent
      \listparindent \itemindent
      \parsep \z@
      }%
   \renewcommand\newblock{\par}}%
}
%
\ExecuteOptions{letterpaper}
\ProcessOptions
%
\newcommand\refname{References}
\newcommand{\ignoretwo}[2]{}
\newcommand{\yearTwoDigits}{\expandafter\ignoretwo\the\year}
\def\@setref#1#2#3{%
  \ifx#1\relax
    \number 0\relax
    \protect\G@refundefinedtrue
    \nfss@text{\reset@font\bfseries ??}%
    \@latex@warning{Reference `#3' on page \thepage \space undefined}%
  \else
    \expandafter#2#1\null
  \fi}
%
%
\lineskip 1pt \normallineskip 1pt
\ifmanuscript
\def\baselinestretch{2}
\else
\def\baselinestretch{1}
\fi
\def\@ixpt{9}
\renewcommand\normalsize{%
  \@setfontsize\normalsize\@ixpt{10pt}
  \abovedisplayskip 6pt plus2pt minus1pt\belowdisplayskip \abovedisplayskip
  \abovedisplayshortskip 6pt plus0pt minus 3pt
  \belowdisplayshortskip 6pt plus0pt minus3pt\let\@listi\@listI} 

\newcommand\small{%
  \@setfontsize\small\@ixpt{11pt}%
  \abovedisplayskip 5pt plus 2pt minus 1pt\belowdisplayskip \abovedisplayskip
  \abovedisplayshortskip 5pt plus0pt minus2pt\belowdisplayshortskip 5pt plus0pt
      minus 2pt
  \def\@listi{\leftmargin\leftmargini \topsep 5pt plus 2pt minus 1pt\parsep 0pt
    plus .7pt 
  \itemsep 1.6pt plus .8pt}}
\newcommand\footnotesize{%
%   \@setfontsize\footnotesize\@viiipt{10pt}
 \@setsize\footnotesize{10pt}\viiipt\@viiipt
  \abovedisplayskip 4pt plus 1pt minus 0pt\belowdisplayskip \abovedisplayskip
  \abovedisplayshortskip 4pt plus 0pt minus 1pt\belowdisplayshortskip 4pt plus
       0pt minus 1pt
  \def\@listi{\leftmargin\leftmargini \topsep 4pt plus 1pt minus
     0pt\parsep 0pt plus .5pt 
     \itemsep 1pt plus .7pt}}

\newcommand\scriptsize{\@setfontsize\scriptsize\@viipt\@viiipt}
\newcommand\tiny{\@setfontsize\tiny\@vpt\@vipt}
\newcommand\large{\@setfontsize\large\@xiipt{14}}
\newcommand\Large{\@setfontsize\Large\@xivpt{18}}
\newcommand\LARGE{\@setfontsize\LARGE\@xviipt{20}}
\newcommand\huge{\@setfontsize\huge\@xxpt{25}}
\newcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}
%
\normalsize 
%
\newdimen\tempdimen  
%
\setlength\trimheight{11in}
\setlength\trimwidth{8.5in}
%
\typeheight52.5pc
\typewidth42pc
\textheight52.5pc
\textwidth42pc
\advance\textheight-3pt
\newdimen\normaltextheight
\setlength\normaltextheight{\textheight}
\oddsidemargin4.5pc
\evensidemargin4.5pc
\topmargin20pt %.25in
\headheight 6pt%
\headsep 29.2pt%
\topskip6pt%
\footskip 100pt 
%
\marginparwidth 0.5in 
\marginparsep .125in 
\columnsep24pt
\columnseprule 0pt 
%
\def\titlefont{\huge\selectfont\centering\mathversion{bold}}
\def\authorfont{\fontfamily{phv}\fontsize{10}{12}\selectfont\rightskip0pt plus1fill} %\mathversion{sfnormal}
\def\rhfont{\fontfamily{phv}\fontsize{9}{10}\selectfont\mathversion{sfnormal}}

\def\sectionfont{\fontfamily{ptm}\fontsize{9}{12}\capsshape\selectfont\raggedright} %\mathversion{rmnormal}
\def\subsectionfont{\fontfamily{ptm}\fontsize{9}{12}\selectfont} %\mathversion{rmnormal}
\def\figcaptionfont{\fontsize{8}{10}\selectfont\mathversion{normal}}%
\def\subcaptionfont{\fontsize{8}{10}\selectfont\mathversion{normal}}%
\def\subcaption#1{{\centering\subcaptionfont#1\par}}
%
\def\tablefont{\fontsize{8}{10}\selectfont}%
\def\tablecaptionfont{\fontsize{9}{11}\selectfont\centering}%
\def\tablenumfont{\fontsize{9}{11}\selectfont}%
\def\tabnotefont{\fontsize{7}{9}\selectfont}
%
\def\encodingdefault{OT1}%
\fontencoding{OT1}%
%
\DeclareFontShape{OMS}{cmsy}{m}{n}{<-> cmsy10 }{}
\DeclareFontShape{OMS}{cmsy}{b}{n}{<-> cmbsy10 }{}
\def\cal{\mathcal}
%
\def\boldmath{\mathversion{bold}}
\def\bm#1{\mathchoice
          {\mbox{\boldmath$\displaystyle#1$}}%
          {\mbox{\boldmath$#1$}}%
          {\mbox{\boldmath$\scriptstyle#1$}}%
          {\mbox{\boldmath$\scriptscriptstyle#1$}}}
%
\footnotesep 7pt
\skip\footins 15pt plus 4pt minus 3pt 
\floatsep 12pt plus 2pt minus 2pt 
\textfloatsep \floatsep 
\intextsep 1pc plus 1pc 
\dblfloatsep 12pt plus 2pt minus 2pt 
\dbltextfloatsep 20pt plus 2pt minus 4pt 
\@fptop 0pt plus 1fil \@fpsep 1pc plus 2fil \@fpbot 0pt plus 1fil 
\@dblfptop 0pt plus 1fil \@dblfpsep 8pt plus 2fil \@dblfpbot 0pt plus 1fil
\marginparpush 6pt 
\parskip 0pt \parindent 0pt \partopsep 0pt  % plus .1pt FBU
\@lowpenalty 51 \@medpenalty 151 \@highpenalty 301 
\@beginparpenalty -\@lowpenalty \@endparpenalty -\@lowpenalty \@itempenalty
-\@lowpenalty 
%
\def\part{\@ucheadtrue
 \@startsection{part}{9}{\z@}{-10pt plus -4pt minus 
 -2pt}{4pt}{\reset@font\normalsize\rmfamily}}
\def\section{\@ucheadtrue
 \@startsection{section}{1}{\z@}{-10pt plus -4pt minus 
 -2pt}{6pt}{\reset@font\fontsize{10}{12}\raggedright\rmfamily\bfseries}}
\def\subsection{\@ucheadfalse
 \@startsection{subsection}{2}{\z@}{-8pt plus -2pt minus 
 -1pt}{6pt}{\reset@font\fontsize{10}{12}\raggedright\rmfamily\bfseries}}
\def\subsubsection{\@ucheadfalse
 \@startsection{subsubsection}{3}{\parindent}{6pt plus 
1pt}{-5pt}{\reset@font\fontsize{9}{10}\itshape}}
\def\paragraph{\@ucheadfalse
 \@startsection{paragraph}{3}{\parindent}{6pt plus 
1pt}{-5pt}{\reset@font\fontsize{10}{12}\itshape}}
%%
\renewcommand{\@seccntformat}[1]{\textup{\csname the#1\endcsname}}
\gdef\@period{.}
\def\@trivlist{\@topsepadd\topsep
\if@noskipsec \gdef\@period{}\leavevmode\gdef\@period{.}\fi
 \ifvmode \advance\@topsepadd\partopsep \else \unskip\par\fi
 \if@inlabel \@noparitemtrue \@noparlisttrue 
 \else \@noparlistfalse \@topsep\@topsepadd \fi
 \advance\@topsep \parskip
 \leftskip\z@\rightskip\@rightskip \parfillskip\@flushglue
 \@setpar{\if@newlist\else{\@@par}\fi} \global\@newlisttrue
\@outerparskip\parskip}
%
\def\@startsection#1#2#3#4#5#6{%
  \if@noskipsec \leavevmode \fi
  \par
  \@tempskipa #4\relax
  \@afterindenttrue
  \ifdim \@tempskipa <\z@
    \@tempskipa -\@tempskipa \@afterindentfalse
  \fi
  \if@nobreak
    \everypar{}%
    \ifnum#2=2
	    \vskip-2pt 
	   \fi
  \else
    \addpenalty\@secpenalty\addvspace\@tempskipa
  \fi
  \@ifstar
    {\@ssect{#3}{#4}{#5}{#6}}%
    {\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}}
%
\def\@sect#1#2#3#4#5#6[#7]#8{%
  \ifnum #2>\c@secnumdepth
    \let\@svsec\@empty
  \else
    \refstepcounter{#1}%
      \if@uchead%
            \protected@edef\@svsec{\@seccntformat{#1}.\quad\relax}%
        \else%
            \protected@edef\@svsec{\@seccntformat{#1}\quad\relax}%
        \fi%
  \fi
  \@tempskipa #5\relax
  \ifdim \@tempskipa>\z@
    \begingroup
      #6{%
        \@hangfrom{\hskip #3\relax\@svsec}%
          \interlinepenalty \@M #8 \@@par}%
    \endgroup
    \csname #1mark\endcsname{#7}%
    \addcontentsline{toc}{#1}{%
      \ifnum #2>\c@secnumdepth \else
        \protect\numberline{\csname the#1\endcsname}%
      \fi
      #7}%
  \else
    \def\@svsechd{%
      #6{\hskip #3\relax
      \@svsec \if@uchead\Makeuppercase{#8}\else#8\fi}%
      \csname #1mark\endcsname{#7}%
      \addcontentsline{toc}{#1}{%
        \ifnum #2>\c@secnumdepth \else
          \protect\numberline{\csname the#1\endcsname}%
        \fi
        #7}}%
  \fi
  \@xsect{#5}}

\def\@xsect#1{\@tempskipa #1\relax
 \ifdim \@tempskipa>\z@
 \par \nobreak
 \vskip \@tempskipa
 \@afterheading
 \else \global\@nobreakfalse \global\@noskipsectrue
 \everypar{\if@noskipsec \global\@noskipsecfalse
 \clubpenalty\@M \hskip -\parindent
 \begingroup \@svsechd\@period \endgroup \unskip
 \hskip -#1
 \else \clubpenalty \@clubpenalty
 \everypar{}\fi}\fi\ignorespaces}
\newif\if@uchead\@ucheadfalse
%
\setcounter{secnumdepth}{3}
\newcounter{secnumbookdepth}
\setcounter{secnumbookdepth}{3}
\newfont{\apbf}{cmbx9}
\def\appendix{\par
 \setcounter{section}{0}
 \setcounter{subsection}{0}
  \section*{APPENDIX}\vskip10pt
        \def\thesection{\Alph{section}}
        \def\theHsection{\Alph{section}}}
%
\labelsep 4pt
\settowidth{\leftmargini}{(9)} \addtolength\leftmargini\labelsep
\settowidth{\leftmarginii}{(b)} \addtolength\leftmarginii\labelsep
\leftmarginiii \leftmarginii
\leftmarginiv \leftmarginii 
\leftmarginv \leftmarginii 
\leftmarginvi \leftmarginii 
\leftmargin\leftmargini
\labelwidth\leftmargini\advance\labelwidth-\labelsep
\def\@listI{\leftmargin\leftmargini \parsep 0pt plus 1pt\topsep 6pt plus 2pt
minus 2pt\itemsep 2pt plus 1pt minus .5pt}
\let\@listi\@listI
\@listi 
\def\@listii{\leftmargin\leftmarginii
 \labelwidth\leftmarginii\advance\labelwidth-\labelsep
 \topsep 0pt plus 1pt 
 \parsep 0pt plus .5pt 
 \itemsep \parsep}
\def\@listiii{\leftmargin\leftmarginiii
 \labelwidth\leftmarginiii\advance\labelwidth-\labelsep
 \topsep 0pt plus 1pt 
 \parsep 0pt plus .5pt 
 \itemsep \parsep}
\def\@listiv{\leftmargin\leftmarginiv
 \labelwidth\leftmarginiv\advance\labelwidth-\labelsep}
\def\@listv{\leftmargin\leftmarginv
 \labelwidth\leftmarginv\advance\labelwidth-\labelsep}
\def\@listvi{\leftmargin\leftmarginvi
 \labelwidth\leftmarginvi\advance\labelwidth-\labelsep}
%
\def\enumerate{\ifnum \@enumdepth >3 \@toodeep\else
 \advance\@enumdepth \@ne 
 \edef\@enumctr{enum\romannumeral\the\@enumdepth}\list
 {\csname label\@enumctr\endcsname}{\usecounter
 {\@enumctr}\def\makelabel##1{##1\hss}}\fi}
\def\longenum{\ifnum \@enumdepth >3 \@toodeep\else
 \advance\@enumdepth \@ne 
 \edef\@enumctr{enum\romannumeral\the\@enumdepth}\list
 {\csname label\@enumctr\endcsname}{\usecounter
 {\@enumctr}\labelwidth\z@}\fi}
\let\endlongenum\endlist
\def\labelenumi{{\rm (}\arabic{enumi}\/{\rm )}} 
\def\theenumi{\arabic{enumi}} 
\def\labelenumii{{\rm (}\alph{enumii}\rm{)}}
\def\theenumii{\alph{enumii}}
\def\p@enumii{\theenumi}
\def\labelenumiii{\roman{enumiii}.}
\def\theenumiii{\roman{enumiii}}
\def\p@enumiii{\theenumi{\rm (}\theenumii{\rm )}}
\def\labelenumiv{\Alph{enumiv}.}
\def\theenumiv{\Alph{enumiv}} 
\renewcommand\theenumiv{\@Alph\c@enumiv}
\def\p@enumiv{\p@enumiii\theenumiii}

\def\p@enumiv{\p@enumiii\theenumiii}

\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}

\def\itemize{\list{---\hskip -\labelsep}{\settowidth
 {\leftmargin}{---}\labelwidth\leftmargin
 \addtolength{\labelwidth}{-\labelsep}}}
\let\enditemize\endlist
\def\longitem{\list{---}{\labelwidth\z@
 \leftmargin\z@ \itemindent\parindent \advance\itemindent\labelsep}}
\let\endlongitem\endlist
\def\verse{\let\\=\@centercr 
 \list{}{\leftmargin 2pc 
 \itemindent -1.5em\listparindent \itemindent 
 \rightmargin\leftmargin\advance\leftmargin 1.5em}\item[]}
\let\endverse\endlist
\def\quotation{\list{}{\leftmargin 2pc \listparindent .5em
 \itemindent\listparindent
 \rightmargin\leftmargin \parsep 0pt plus 1pt}\item[]}
\let\endquotation=\endlist
\def\quote{\list{}{\leftmargin 2pc \rightmargin\leftmargin}\item[]}
\let\endquote=\endlist

%
\newenvironment{unnumlist}{%
    \list{}{%
    \listparindent\parindent
    \itemindent-1em
    \leftmargin1em
    \parsep0pt
    \itemsep2pt
    \partopsep0pt}
    \def\makelable##1{##1}%
}{\endlist}%
%
\def\description{\list{}{\listparindent\parindent\labelwidth\z@
 \leftmargin\z@ \itemindent\parindent\advance\itemindent\labelsep
 \def\makelabel##1{\it ##1}}}
\let\enddescription\endlist
%
\def\describe#1{\list{}{\listparindent\parindent\settowidth{\labelwidth}{#1}\leftmargin
 \labelwidth\addtolength\leftmargin\labelsep\def\makelabel##1{##1\hfil}}}
\let\enddescribe\endlist
%
\def\program{\ifx\@currsize\normalsize\small \else \rm \fi\tabbing}
\let\endprogram\endtabbing
%
\newtheorem{theorem}{Theorem}
\newtheorem{strategy}{Strategy}
\newtheorem{property}{Property}
\newtheorem{proposition}{Proposition}
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{exam}{Example}
\newenvironment{example}{%
\italicenvfalse
\begin{exam}}{\end{exam}\italicenvtrue}
%
\newtheorem{defi}[theorem]{Definition}
\newenvironment{definition}{%
\italicenvfalse
\begin{defi}}{\end{defi}\italicenvtrue}
%
\def\@begintheorem#1#2{\trivlist \item[\hskip 10pt\hskip 
   \labelsep{\sc{#1}\hskip 5pt\relax #2.}] \itshape}
%
\def\@opargbegintheorem#1#2#3{\trivlist
     \item[\hskip 10pt \hskip 
\labelsep{\sc{#1}\savebox\@tempboxa{\sc{#3}}\ifdim 
        \wd\@tempboxa>\z@ \hskip 5pt\relax \sc{#2}  \box\@tempboxa\fi.}]
\itshape}
%
\newif\if@qeded\global\@qededfalse
\def\proof{\global\@qededfalse\@ifnextchar[{\@xproof}{\@proof}}
\def\endproof{\if@qeded\else\qed\fi\endtrivlist}
\def\qed{\unskip\kern 10pt{\unitlength1pt\linethickness{.4pt}\framebox(5,5){}}
\global\@qededtrue}
\def\@proof{\trivlist \item[\hskip 10pt\hskip 
 \labelsep{\sc Proof.}]\ignorespaces}
\def\@xproof[#1]{\trivlist \item[\hskip 10pt\hskip 
 \labelsep{\sc Proof #1.}]\ignorespaces}
%
\def\newdef#1#2{\expandafter\@ifdefinable\csname #1\endcsname
{\@definecounter{#1}\expandafter\xdef\csname
the#1\endcsname{\@thmcounter{#1}}\global
 \@namedef{#1}{\@defthm{#1}{#2}}\global
 \@namedef{end#1}{\@endtheorem}}}
\def\@defthm#1#2{\refstepcounter
 {#1}\@ifnextchar[{\@ydefthm{#1}{#2}}{\@xdefthm{#1}{#2}}}
\def\@xdefthm#1#2{\@begindef{#2}{\csname the#1\endcsname}\ignorespaces}
\def\@ydefthm#1#2[#3]{\trivlist \item[\hskip 10pt\hskip 
 \labelsep{\it #2\savebox\@tempboxa{#3}\ifdim 
 \wd\@tempboxa>\z@ \ \box\@tempboxa\fi.}]\ignorespaces}
\def\@begindef#1#2{\trivlist \item[\hskip 10pt\hskip 
 \labelsep{\it #1\ \rm #2.}]}
%
\def\theequation{\arabic{equation}}
%
\def\titlepage{\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
 \else \newpage \fi \thispagestyle{empty}\c@page\z@}
\def\endtitlepage{\if@restonecol\twocolumn \else \newpage \fi}
%
\arraycolsep 2.5pt \tabcolsep 6pt \arrayrulewidth .4pt \doublerulesep 2pt 
\tabbingsep \labelsep 
%
\skip\@mpfootins = \skip\footins
\fboxsep = 3pt \fboxrule = .4pt 
%
\newcounter{part}
\newcounter{section}
\newcounter{subsection}[section]
\newcounter{subsubsection}[subsection]
\newcounter{paragraph}[subsubsection]
%
\def\thepart{\Roman{part}} 
\def\thesection {\arabic{section}}
\def\thesubsection {\thesection.\arabic{subsection}}
\def\thesubsubsection {\itshape\thesubsection.\arabic{subsubsection}}
\def\theparagraph {\thesubsubsection.\arabic{paragraph}}

\def\@pnumwidth{1.55em}
\def\@tocrmarg {2.55em}
\def\@dotsep{4.5}
\setcounter{tocdepth}{3}

\def\tableofcontents{\section*{Contents\@mkboth{CONTENTS}{CONTENTS}}
 \@starttoc{toc}}
\def\l@part#1#2{\addpenalty{\@secpenalty}
 \addvspace{2.25em plus 1pt} \begingroup
 \@tempdima 3em \parindent \z@ \rightskip \@pnumwidth \parfillskip
-\@pnumwidth 
 {\large \bf \leavevmode #1\hfil \hbox to\@pnumwidth{\hss #2}}\par
 \nobreak \endgroup}
\def\l@section#1#2{\addpenalty{\@secpenalty} \addvspace{1.0em plus 1pt}
\@tempdima 1.5em \begingroup
 \parindent \z@ \rightskip \@pnumwidth 
 \parfillskip -\@pnumwidth 
 \bf \leavevmode #1\hfil \hbox to\@pnumwidth{\hss #2}\par
 \endgroup}
\def\l@subsection{\@dottedtocline{2}{1.5em}{2.3em}}
\def\l@subsubsection{\@dottedtocline{3}{3.8em}{3.2em}}
\def\listoffigures{\section*{List of Figures\@mkboth
 {LIST OF FIGURES}{LIST OF FIGURES}}\@starttoc{lof}}
\def\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}
\def\listoftables{\section*{List of Tables\@mkboth
 {LIST OF TABLES}{LIST OF TABLES}}\@starttoc{lot}}
\let\l@table\l@figure
%
\newif\if@restonecol
\def\theindex{\@restonecoltrue\if@twocolumn\@restonecolfalse\fi
\columnseprule \z@
\columnsep 35pt\twocolumn[\section*{Index}]
 \@mkboth{INDEX}{INDEX}\thispagestyle{plain}\parindent\z@
 \parskip\z@ plus .3pt\relax\let\item\@idxitem}
\def\@idxitem{\par\hangindent 40pt}
\def\subitem{\par\hangindent 40pt \hspace*{20pt}}
\def\subsubitem{\par\hangindent 40pt \hspace*{30pt}}
\def\endtheindex{\if@restonecol\onecolumn\else\clearpage\fi}
\def\indexspace{\par \vskip 10pt plus 5pt minus 3pt\relax}
%
\def\footnoterule{\kern-3\p@
 \hrule \@height 0.2\p@ \@width 47\p@ 
  \kern 2.6\p@
} 

\long\def\@makefntext#1{\parindent 1em\noindent 
 $^{\@thefnmark}$#1}
%
\setcounter{topnumber}{3}
\def\topfraction{.99}
\setcounter{bottomnumber}{1}
\def\bottomfraction{.5}
\setcounter{totalnumber}{3}
\def\textfraction{.01}
\def\floatpagefraction{.85}
\setcounter{dbltopnumber}{2}
\def\dbltopfraction{.95}
\def\dblfloatpagefraction{.96}
%
\long\def\@makecaption#1#2{\vskip 1pc \setbox\@tempboxa\hbox{#1.\hskip
1em\relax #2}
 \ifdim \wd\@tempboxa >\hsize #1. #2\par \else \hbox
to\hsize{\hfil\box\@tempboxa\hfil} 
 \fi}

\def\nocaption{\refstepcounter\@captype \par
 \vskip 1pc \hbox to\hsize{\hfil \footnotesize Figure \thefigure
 \hfil}}
%
\newcounter{figure}
\def\thefigure{\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{Fig.\ \thefigure}%
\def\figure{\let\normalsize\footnotesize\normalsize\@float{figure}}
\let\endfigure\end@float
\@namedef{figure*}{\@dblfloat{figure}}
\@namedef{endfigure*}{\end@dblfloat}
%
\newcounter{table}
\def\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\newlength\belowcaptionskip
\setlength\belowcaptionskip{1\p@}
%
\def\FigName{figure}%
\long\def\@caption#1[#2]#3{\par\begingroup
    \@parboxrestore
    \normalsize \bf \centering
    \@makecaption{\csname fnum@#1\endcsname}{\ignorespaces #3}\par
  \endgroup}
%
%
\newbox\tbbox
\long\def\@makecaption#1#2{%
  \ifx\FigName\@captype
    \vskip 7.3pt
    \setbox\@tempboxa\hbox{\figcaptionfont{#1}.\hskip7.3pt\relax #2\par}%
    \ifdim \wd\@tempboxa >\hsize
        \figcaptionfont{#1}.\hskip7.3pt\relax #2\par
    \else
        \centerline{\box\@tempboxa}%
    \fi
  \else%
    \setbox\tbbox=\vbox{\hsize\tempdimen{\tablenumfont #1}\ {\tablecaptionfont #2\par}}%
    \setbox\@tempboxa\hbox{\hsize\tempdimen{\tablenumfont #1}\ {\tablecaptionfont #2\par}\vphantom{jgq}}%
    \ifdim \wd\@tempboxa >\tempdimen
        \centerline{\box\tbbox}%        
    \else
        \centerline{\box\@tempboxa}%
    \fi
    \vskip\belowcaptionskip
  \fi}
%
\def\fnum@table{Table~\thetable.\ }
\def\table{\let\normalsize\footnotesize \normalsize\@float{table}}
\let\endtable\end@float
\@namedef{table*}{\@dblfloat{table}}
\@namedef{endtable*}{\end@dblfloat}
\def\ijcatable#1{\@narrowfig #1\relax
 \let\caption\@atcap \let\nocaption\@atnocap
 \def\@tmpnf{}\@ifnextchar[{\@xntab}{\@ntab}}
\def\endijcatable{\hbox to \textwidth{\hfil
\vbox{\hsize \@narrowfig
\box\@nfcapbox
{\baselineskip 4pt \hbox{\vrule height .4pt width \hsize}}
\vskip -1pt
\box\@nfigbox\vskip -1pt
{\baselineskip 4pt \hbox{\vrule height .4pt width \hsize}}}\hfil}
\end@float}
\def\@xntab[#1]{\def\@tmpnf{[#1]}\@ntab}
\def\@ntab{\expandafter\table\@tmpnf
 \setbox\@nfigbox\vbox\bgroup
 \hsize \@narrowfig \@parboxrestore}
\def\@atmakecap #1#2{\setbox\@tempboxa\hbox{#1.\hskip 1em\relax #2}
 \ifdim \wd\@tempboxa >\hsize \sloppy #1.\hskip 1em\relax #2 \par \else \hbox
to\hsize{\hfil\box\@tempboxa\hfil}
 \fi}
\def\@atcap{\par\egroup\refstepcounter\@captype 
 \@dblarg{\@atcapx\@captype}}
\long\def\@atcapx#1[#2]#3{\setbox\@nfcapbox\vbox {\hsize \wd\@nfigbox
 \@parboxrestore
 \@atmakecap{\csname fnum@#1\endcsname}{\ignorespaces #3}\par}}
\def\@atnocap{\egroup \refstepcounter\@captype 
 \setbox\@nfcapbox\vbox {\hsize \wd\@nfigbox
 \hbox to\hsize{\hfil \footnotesize Table \thetable\hfil}}}
%
\newdimen\tabledim
%
\long\def\tbl#1#2{%
 \setbox\tempbox\hbox{\tablefont #2}%
 \tabledim\hsize\advance\tabledim by -\wd\tempbox
 \tempdimen\wd\tempbox
	\global\divide\tabledim\tw@
 \caption{#1}
	\centerline{\box\tempbox}
  }%
% 
\newenvironment{tabnote}{%
\par%\addvspace{-1pt}
\tabnotefont
\@ifnextchar[{\@tabnote}{\@tabnote[]}}{%
\par}
\def\@tabnote[#1]{\def\@Tempa{#1}\leftskip\tabledim\rightskip\leftskip\ifx\@Tempa\@empty\else{\it #1:}\ \fi\ignorespaces}
%
\def\tabnoteentry#1#2{\parindent0pt\par\@hangfrom{#1}{#2}}
\def\Note#1#2{\parindent0pt\par\hangindent3.7pt{\it #1}\ #2}
%

\def\Hline{%
  \noalign{\ifnum0=`}\fi\hrule \@height .5pt \futurelet
   \@tempa\@xhline}
%
\def\narrowfig#1{\@narrowfig #1\relax
 \let\caption\@nfcap \let\nocaption\@nfnocap
 \def\@tmpnf{}\@ifnextchar[{\@xnfig}{\@nfig}}
\def\endnarrowfig{\hbox to \textwidth{\if@nfeven
 \box\@nfcapbox\hfil\box\@nfigbox
 \else \box\@nfigbox\hfil\box\@nfcapbox\fi}\end@float}
\def\@xnfig[#1]{\def\@tmpnf{[#1]}\@nfig}
\def\@nfig{\expandafter\figure\@tmpnf
 \setbox\@nfigbox\vbox\bgroup
 \hsize \@narrowfig \@parboxrestore}
\def\@nfmakecap #1#2{\setbox\@tempboxa\hbox{#1.\hskip 1em\relax #2}
 \ifdim \wd\@tempboxa >\hsize \sloppy #1.\hskip 1em\relax #2 \par \else \hbox
to\hsize{\if@nfeven\else\hfil\fi\box\@tempboxa\if@nfeven\hfil\fi}
 \fi}
\def\@nfcap{\par\egroup\refstepcounter\@captype 
 \@dblarg{\@nfcapx\@captype}}
\long\def\@nfcapx#1[#2]#3{\@seteven
 \setbox\@nfcapbox\vbox to \ht\@nfigbox
 {\hsize \textwidth \advance\hsize -2pc \advance\hsize -\wd\@nfigbox
 \@parboxrestore
 \vfil 
 \@nfmakecap{\csname fnum@#1\endcsname}{\ignorespaces #3}\par
 \vfil}}
\def\@nfnocap{\egroup \refstepcounter\@captype \@seteven
 \setbox\@nfcapbox\vbox to \ht\@nfigbox
 {\hsize \textwidth \advance\hsize -2pc \advance\hsize -\wd\@nfigbox
 \@parboxrestore
 \vfil 
 \hbox to\hsize{\if@nfeven\else\hfil\fi 
 \footnotesize Figure \thefigure
 \if@nfeven\hfil\fi}
 \vfil}}
\def\@seteven{\@nfeventrue 
   \@ifundefined{r@@nf\thefigure}{}{%
      \edef\@tmpnf{\csname r@@nf\thefigure\endcsname}%
      \edef\@tmpnf{\expandafter\@getpagenum\@tmpnf}%
      \ifodd\@tmpnf\relax\@nfevenfalse\fi}%
\label{@nf\thefigure}\edef\@tmpnfx{\if@nfeven e\else o\fi}
\edef\@tmpnf{\write\@unused {\noexpand\ifodd \noexpand\c@page 
 \noexpand\if \@tmpnfx e\noexpand\@nfmsg{\thefigure} \noexpand\fi
 \noexpand\else 
 \noexpand\if \@tmpnfx o\noexpand\@nfmsg{\thefigure}\noexpand\fi 
 \noexpand\fi }}\@tmpnf}
\def\@nfmsg#1{Bad narrowfig: Figure #1 on page \thepage}

\newdimen\@narrowfig
\newbox\@nfigbox
\newbox\@nfcapbox
\newif\if@nfeven 


\def\maketitle{%
   \thispagestyle{titlepage}%
     \newpage
     \global\@topnum\z@
     \twocolumn[\@maketitle]%
    \let\maketitle\relax
    \global\let\@sponsors\@empty
}
%
\def\@maketitle{\newpage \thispagestyle{titlepage}\par
 \begingroup \lineskip = \z@\null
 \begin{picture}(5,5)
	\includegraphics[width=1in]{logojuliacon.pdf}
 \end{picture}
 \vskip -7pt\relax %-18.5pt
 \parindent\z@ \LARGE {\centering \hyphenpenalty\@M
  {\titlefont \@title} \par
 \global\firstfoot %aiellom
 \global\runningfoot %aiellom
}
\label{@firstpg}
{ 
\begin{center}%
    \vskip 0.1em%
    {\large
     \lineskip .75em%
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
      \vskip 1.5em%
  \end{center}\par
  \@thanks
}
 \vskip 23pt\relax
 \endgroup
 }
\newbox\@abstract
\newbox\@terms
\newbox\@keywords

%
\newenvironment{abstract}
	       {\section*{ABSTRACT}\par\fontsize{10}{12}\indent\ignorespaces}
	       {
		   { \ifvoid\@terms\else\box\@terms\fi
 \@keywords \@juliaconformat\empty}\vskip6pt}
%
\def\terms#1{\setbox\@terms=\vbox{\hsize20pc%
 \footnotesize%
\parindent 0pt \noindent 
 { \section*{General Terms}} \ignorespaces #1}}
\def\keywords#1{\gdef\@keywords{\hsize20pc%
\parindent 0pt\noindent\ignorespaces% 
 { \section*{Keywords} } \ignorespaces #1}
}

\def\category#1#2#3{\@ifnextchar
 [{\@category{#1}{#2}{#3}}{\@xcategory{#1}{#2}{#3}}}
\def\@category#1#2#3[#4]{\edef\@tempa{\ifx \@categories\@empty 
 \else ; \fi}{\def\protect{\noexpand\protect
 \noexpand}\def\and{\noexpand\and}\xdef\@categories{\@categories\@tempa #1
[{\bf #2}]: 
 #3\kern\z@---\hskip\z@{\it #4}}}}
\def\@xcategory#1#2#3{\edef\@tempa{\ifx \@categories\@empty \else ;
\fi}{\def\protect{\noexpand\protect\noexpand}\def\and{\noexpand
 \and}\xdef\@categories{\@categories\@tempa #1 [{\bf #2}]: #3}}}
\def\@categories{}

\newenvironment{ackslike}[1]
  {\par \footnotesize 
   \@ucheadfalse
   \@startsection{subsection}{2}{\z@}{-16pt plus -2pt minus -1pt}{2pt}{\sf}*
   {\uppercase{#1}}\par\normalsize
        }
  {\par}
\newenvironment{acks}{\begin{ackslike}{ \normalsize\rm\bf Acknowledgments}}{\end{ackslike}}
%

\newcommand\headingtable{%
  \begin{tabular}[b]{l} {\@journalName}\end{tabular}}
\markright{\protect\headingtable}
\mark{{}{}} 
\def\bull{{\fontsize{7}{7}\selectfont\raise1.6pt\hbox{$\bullet$}}}
\def\ps@myheadings{\let\@mkboth\@gobbletwo
\def\@oddhead{ \fontsize{9}{12} \rm {{\itshape\headingtable}\hfill \@volume(\@issue), \@year}}
\def\@oddfoot{\fontsize{9}{12}\@runningfoot}
\def\@evenhead{ \fontsize{9}{12} \rm {\itshape\headingtable}\hfill \@volume(\@issue), \@year}
\def\@evenfoot{\fontsize{9}{12}\@runningfoot}
\def\sectionmark##1{}\def\subsectionmark##1{}}
%
\def\@runningfoot{}
\def\runningfoot{\def\@runningfoot{ \fontsize{9}{12} \thepage}}
\def\@firstfoot{}
\def\firstfoot{\def\@firstfoot{\fontsize{9}{12} \thepage}}
\def\ps@titlepage{\let\@mkboth\@gobbletwo
\def\@oddhead{}\def\@oddfoot{\fontsize{9}{12}\@firstfoot}\def\@evenhead{}\def\@evenfoot{\fontsize{9}{12}\@firstfoot}}
%
\def\today{\ifcase\month\or
 January\or February\or March\or April\or May\or June\or
 July\or August\or September\or October\or November\or December\fi
 \space\number\day, \number\year}
\def\@marrayclassiv{\@addtopreamble{$\displaystyle \@nextchar$}}
\def\@marrayclassz{\ifcase \@lastchclass \@acolampacol \or \@ampacol \or
 \or \or \@addamp \or
 \@acolampacol \or \@firstampfalse \@acol \fi
\edef\@preamble{\@preamble
 \ifcase \@chnum 
 \hfil$\relax\displaystyle\@sharp$\hfil \or $\relax\displaystyle\@sharp$\hfil 
 \or \hfil$\relax\displaystyle\@sharp$\fi}}
\def\marray{\arraycolsep 2.5pt\let\@acol\@arrayacol \let\@classz\@marrayclassz
 \let\@classiv\@marrayclassiv \let\\\@arraycr\def\@halignto{}\@tabarray}
\def\endmarray{\crcr\egroup\egroup}
%
\ps@myheadings \pagenumbering{arabic} \onecolumn 
%
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
%
\def\newdef#1{\@ifnextchar[{\@xnewdef{#1}}{\@ynewdef{#1}}}
\def\@xnewdef#1[#2]#3{\newtheorem{italic@#1}[#2]{{\em #3}}\@newdef{#1}}
\def\@ynewdef#1#2{\@ifnextchar[{\@xynewdef{#1}{#2}}{\@yynewdef{#1}{#2}}}
\def\@xynewdef#1#2[#3]{\newtheorem{italic@#1}{{\em #2}}[#3]\@newdef{#1}}
\def\@yynewdef#1#2{\newtheorem{italic@#1}{{\em #2}}\@newdef{#1}}
\def\@newdef#1{\newenvironment{#1}{\@ifnextchar[{\@xstartdef{#1}}{\@ystartdef{#1}}}{\end{italic@#1}}}
\def\@xstartdef#1[#2]{\begin{italic@#1}[{\em #2}]\rm}
\def\@ystartdef#1{\begin{italic@#1}\rm}
%
%\def\@oddfoot{\hbox{}\hfill\@runningfoot \thepage}
%\def\@evenfoot{\@runningfoot\hfill\hbox{} \thepage }
%\def\firstfootsize{\@setsize\firstfootsize{9pt}\vipt\@vipt}
\def\ps@titlepage{\let\@mkboth\@gobbletwo
\def\@oddhead{\fontsize{9}{12} \rm {\hskip 19pt\itshape\headingtable}}\def\@oddfoot{\hbox{}\hfill\fontsize{9}{12}\@firstfoot}%
\def\@evenhead{}\def\@evenfoot{\firstfootsize\@firstfoot\hfill\hbox{}}}
%
\def\@listI{\leftmargin\leftmargini
 \labelwidth\leftmargini\advance\labelwidth-\labelsep
 \parsep 0pt plus 1pt
 \topsep 6pt plus 2pt minus 2pt
 \itemsep 2pt plus 1pt minus .5pt}
\let\@listi\@listI
\@listi 
%
\def\longenum{\ifnum \@enumdepth >3 \@toodeep\else
 \advance\@enumdepth \@ne 
 \edef\@enumctr{enum\romannumeral\the\@enumdepth}\list
 {\csname label\@enumctr\endcsname}{\usecounter
 {\@enumctr}\labelwidth\z@\leftmargin\z@
 \itemindent\parindent \advance\itemindent\labelsep}\fi}
%
\def\ack{ \par \footnotesize 
\@ucheadfalse
\@startsection{subsection}{2}{\z@}{-16pt plus -2pt minus 
 -1pt}{2pt}{\sf}*{ACKNOWLEDGMENT}\par\normalsize
}
\def\endack{\par}

% provide both spellings of Acknowledgment(s)
\let\acknowledgments\acks
\let\endacknowledgments\endacks
\let\acknowledgment\ack
\let\endacknowledgment\endack
%
\newcommand{\bibemph}[1]{{\em#1}}
\newcommand{\bibemphic}[1]{{\em#1\/}}
\newcommand{\bibsc}[1]{{\sc#1}}

\newcommand\bibyear[2]{%
    \unskip{\hskip8pt}\ignorespaces#1\unskip
    \if..#2{\hskip6pt}\else {\hskip8pt}#2 \fi
}
%
\let\l@table\l@figure
\newdimen\bibindent
\setlength\bibindent{1.5em}
\newenvironment{thebibliography}[1]
     {\section{\refname}%%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty

%
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
%
\def\@juliaconformat{}
\def\juliaconformat#1{\gdef\@juliaconformat{\noindent{\bf JuliaCon Reference Format:}\\[2pt] #1\par}}
%
\def\received#1#2{%
    \par%
    \tiny
    \addvspace{12\p@}%
    \parindent\z@%
\small\scriptsize{Received\ #1;\   accepted\ #2}%
\par%
}

%
\sloppy
\clubpenalty10000
\widowpenalty10000%
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
%
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty

\voffset-5pc
\hoffset-6.03pc

\usepackage{times}
%% \usepackage[mtbold]{mathtime}
\usepackage{bm}
\usepackage{graphicx}% Include figure files
\usepackage{hyperref}
%%\usepackage{microtype}
\renewcommand{\ttdefault}{cmtt}

\usepackage{jlcode}

\usepackage{authblk}

\endinput

% end of juliacon.cls