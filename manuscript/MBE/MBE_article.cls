%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                          %%%%%%
%%%%%              * CIPL LATEX CLASS FILE FOR JOURNALS *      %%%%%%
%%%%%                                                          %%%%%%
%%%%%   Custumer: OUP                                          %%%%%%
%%%%%   Journal:  Molbev (Molbev)                              %%%%%%
%%%%%                                                          %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{MBE}
              [2007/03/23 v1.4b
 Standard LaTeX document class]
\newif\if@restonecol
\newif\if@mstype\@mstypefalse
\newif\if@fpage\@fpagefalse
\newif\if@bookreview\@bookreviewfalse
\newif\if@titlepage\@titlepagefalse
\newif\if@openright
\newif\if@mainmatter\@mainmattertrue
\DeclareOption{a4paper}
   {\setlength\paperheight {11in}%
    \setlength\paperwidth  {8.25in}}
\DeclareOption{a5paper}
   {\setlength\paperheight {210mm}%
    \setlength\paperwidth  {148mm}}
\DeclareOption{b5paper}
   {\setlength\paperheight {250mm}%
    \setlength\paperwidth  {176mm}}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{legalpaper}
   {\setlength\paperheight {14in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{executivepaper}
   {\setlength\paperheight {10.5in}%
    \setlength\paperwidth  {7.25in}}
\DeclareOption{custompaper}
   {\setlength\paperheight {276mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{landscape}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
%B%\DeclareOption{mynatbib}{\@mynatbibtrue}
%B%\DeclareOption{nomynatbib}{\@mynatbibfalse}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
\newif\if@cropmarkson \@cropmarksontrue
%\newif\if@mynatbib \@mynatbibtrue
\newif\if@theoremsty \@theoremstytrue
\newif\if@surnamecolor \@surnamecolorfalse
\newif\if@CE \@CEfalse
\newif\if@MSC \@MSCfalse
\DeclareOption{b4center}{\PassOptionsToPackage{b4center}{crop}}
\DeclareOption{a4center}{\PassOptionsToPackage{a4center}{crop}}
\DeclareOption{surname}{\@surnamecolortrue}
\DeclareOption{nosurname}{\@surnamecolorfalse}
\DeclareOption{FPAGE}{\@fpagetrue}
\DeclareOption{CE}{\@CEtrue}
\DeclareOption{noCE}{\@CEfalse}
\DeclareOption{PV}{\@mstypetrue}
\DeclareOption{BR}{\@bookreviewtrue}
\DeclareOption{MSC}{\@MSCtrue}
\DeclareOption{noMSC}{\@MSCfalse}
\DeclareOption{reviews}{\PassOptionsToPackage{cam}{crop}\global\@cropmarksontrue}
\DeclareOption{crop}{\PassOptionsToPackage{cam}{crop}\global\@cropmarksontrue}
\DeclareOption{nocrop}{\PassOptionsToPackage{off}{crop}\global\@cropmarksonfalse}
\DeclareOption{grid}{\setlength\overfullrule{0pt}\AtBeginDocument{\DRAWGRID}}
\DeclareOption{proof}{\gdef\artdir{//cephaserver1/Artwork/Journals/Inprocess/OUP/}}
\DeclareOption{revises}{\gdef\artdir{//cephaserver1/Artwork/Journals/Inprocess/OUP/}}
\DeclareOption{oneside}{\@twosidefalse \@mparswitchfalse}
\DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}
\DeclareOption{draft}{\setlength\overfullrule{5pt}}
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\DeclareOption{webpdf}{\gdef\@webpdf{webpdf}}
\DeclareOption{titlepage}{\@titlepagetrue}
\DeclareOption{notitlepage}{\@titlepagefalse}
\DeclareOption{openright}{\@openrighttrue}
\DeclareOption{openany}{\@openrightfalse}
\DeclareOption{onecolumn}{\@twocolumnfalse}
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{notheoremsty}{\@theoremstyfalse}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\newif\if@openbib \@openbibfalse
\DeclareOption{openbib}{%
\@openbibtrue%
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
\ExecuteOptions{custompaper,twoside,crop,b4center}
\ProcessOptions
\RequirePackage{crop,inputenc,amsmath,amssymb,latexsym,
                float,epsfig,wrapfig,graphicx,stfloats,
                multicol,rotating,multirow,dcolumn,
                boxedminipage,makeidx,soul,url,xspace,color}
%03-12-13%\RequirePackage[LY1]{fontenc}
%03-12-13%\RequirePackage[LY1,mtbold]{MBE-fonts}
%
%\usepackage[mathlines,switch,modulo,pagewise]{lineno}%\pagewiselinenumbers[switch]
%
%\modulolinenumbers[5]
%%%%%%%%%%%%%%%%%%%%%%%% *** SIZES *** %%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%
\renewcommand\normalsize{%
   \@setfontsize\normalsize{10.5\p@}{22\p@}
   %\@setfontsize\normalsize\@xpt\@xiipt
   \abovedisplayskip 5.5\p@ \@plus2\p@ \@minus0\p@
   \abovedisplayshortskip \z@ \@plus3\p@
   \belowdisplayshortskip 3\p@ \@plus3\p@ \@minus0\p@
   \belowdisplayskip \abovedisplayskip
   \let\@listi\@listI}
\normalsize
\newcommand\small{%
   \@setfontsize\small\@ixpt{11}%
   \abovedisplayskip 9\p@ \@plus3\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus2\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip}

\newcommand\footnotesize{%
   \@setfontsize\footnotesize\@viiipt{10}%
   \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus\p@
   \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 3\p@ \@plus\p@ \@minus\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip}
%
\newcommand\scriptsize{\@setfontsize\scriptsize\@viipt\@viiipt}
\newcommand\tiny{\@setfontsize\tiny\@vpt\@vipt}
\newcommand\large{\@setfontsize\large\@xiipt{14}}
\newcommand\Large{\@setfontsize\Large\@xivpt{16}}
\newcommand\LARGE{\@setfontsize\LARGE\@xviipt{22}}
\newcommand\huge{\@setfontsize\huge\@xxpt{25}}
\newcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}
%
%\DeclareMathSizes{text-size}{mathtext-size}{script-size}{scriptscript-size}
%
\DeclareMathSizes{10}{10}{8}{6}
\DeclareMathSizes{9}{9}{7.6}{5}
\DeclareMathSizes{8}{8}{6}{5}

%sen%%\DeclareMathSizes{5}     {5}   {5}  {5}
%sen%\DeclareMathSizes{6}     {6}   {5}  {5}
%sen%\DeclareMathSizes{7}     {7}   {5}  {5}
%sen%\DeclareMathSizes{8}     {8}   {6}  {5}
%sen%\DeclareMathSizes{9}     {9}   {6.5}  {5}
%sen%\DeclareMathSizes{10}   {10}   {7.5}  {5}
%sen%\DeclareMathSizes{12}   {12}   {9}  {7}
%%%%%%%%%%%%%%%%%%%% ** PENALTIES ** %%%%%%%%%%%%%%%%%%%
%
% Page break penalties
\@lowpenalty   31
\@medpenalty  151
\@highpenalty 301
% Disallow widows and orphans
\clubpenalty 10000
\widowpenalty 10000
% Disable page breaks before equations, allow pagebreaks after
% equations and discourage widow lines before equations.
\displaywidowpenalty 100
\predisplaypenalty   10000
\postdisplaypenalty  100
% Allow breaking the page in the middle of a paragraph
\interlinepenalty 0
% Disallow breaking the page after a hyphenated line
\brokenpenalty 10000
% Hyphenation; don't split words into less than three characters
\lefthyphenmin=3
\righthyphenmin=3
\hyphenation{figure figures Figure table tables Table Section section Chapter chapter}
%
% Prevent Hyphenation
\hyphenpenalty=2000
\tolerance=1000
%
\frenchspacing
\vbadness=9999
\doublehyphendemerits=10000
\doublehyphendemerits 640000   % corresponds to badness 800
\finalhyphendemerits  1000000  % corresponds to badness 1000
%
%%%%%%%%%%%%%%%%%%%% ** PAGE LAYOUT ** %%%%%%%%%%%%%%%%%%%
%
% Set page layout dimensions
\setlength\topmargin{13pt}          % head margin
\addtolength\topmargin{-1in}       % subtract out the 1 inch driver margin
%
\setlength\headheight{7\p@}         % height of running head
\setlength\topskip{10\p@}           % height of first line of text
\setlength\headsep{37\p@}         % space below running head --
\addtolength\headsep{-\topskip}     % base to base with first line of text
%
\setlength\textwidth{41.4pc}           % text measure excluding margins
%%\setlength\textheight{58\baselineskip} % # lines of text%%changed by Bala 06-11-2013
%\setlength\textheight{30.5\baselineskip} % # lines of text%%changed by Bala 06-11-2013
\setlength\textheight{31.6\baselineskip} % # lines of text
%\addtolength\textheight{\topskip}
\addtolength\textheight{2pt}
%
% Side Margins
\setlength\oddsidemargin{50.6pt}
\addtolength\oddsidemargin{-1in}   % subtract out the 1 inch driver margin
%
\setlength\evensidemargin{\paperwidth}
\addtolength\evensidemargin{-\textwidth}
\addtolength\evensidemargin{-\oddsidemargin}
\addtolength\evensidemargin{-2in}
%
\setlength\columnsep{24\p@}         % space between columns for double-column text
\setlength\columnseprule{0\p@}      % width of rule between two columns
%
\setlength\lineskip{1\p@}
\setlength\normallineskip{10\p@}
\renewcommand\baselinestretch{1}
\setlength\parskip{0\p@}
\setlength\parindent{1em}
\parindent=10pt
%
\setlength\footskip{20\p@}
\setlength\maxdepth{0pt}%
%
\setlength\marginparwidth{2.5pc}
\setlength\marginparsep{.5pc}
\def \@marginparreset {%
        \reset@font
        \scriptsize
        \@setminipage}
%
\setlength\smallskipamount{3\p@ \@plus 1\p@ \@minus 1\p@}
\setlength\medskipamount{6\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\bigskipamount{12\p@ \@plus 4\p@ \@minus 4\p@}
%
%%%%%%%%%%%%%%%%%%%% ** TITLE DEFS ** %%%%%%%%%%%%%%%%%%%
%
\newcommand\contentsname{Contents}
\newcommand\listfigurename{List of Figures}
\newcommand\listtablename{List of Tables}
\newcommand\refname{{References}}
\newcommand\chaptername{Chapter}%
\newcommand\indexname{Index}
\newcommand\figurename{{F\textsc{IG.}}}
\newcommand\tablename{{Table}}
\newcommand\partname{Part}
\newcommand\appendixname{Appendix}
\newcommand\abstractname{{\abstractheadfont Abstract}}
\newcommand\ackname{Acknowledgements}
\newcommand\keywordname{Key{\hskip3.5pt}words}
\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
%
%%%%%%%%%%%%%%%%%%%% ** OLD COMMANDS ** %%%%%%%%%%%%%%%%%%%
%
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\sfit}{\normalfont\sffamily\itshape}{\mathsf}
\DeclareOldFontCommand{\sfb}{\normalfont\sffamily\bfseries}{\mathsf}
\DeclareOldFontCommand{\sfbi}{\normalfont\sffamily\bfseries\itshape}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\tti}{\normalfont\ttfamily\itshape}{\mathtt}
\DeclareOldFontCommand{\ttb}{\normalfont\ttfamily\bfseries}{\mathtt}
\DeclareOldFontCommand{\ttbi}{\normalfont\ttfamily\bfseries\itshape}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\bi}{\normalfont\bfseries\itshape}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
%
%%%%%%%%%%%%%%%%%%%% ** PAGE STYLES ** %%%%%%%%%%%%%%%%%%%
%
\def\volname#1{\gdef\@volname{#1}}
\volname{VOL.}
%
\if@twoside
  \def\ps@headings{%
      %\let\@oddfoot\@empty
      %\let\@evenfoot\@empty
      \def\@evenhead{\GRIDRULES\llap{\hrule height .5pt width \textwidth}{\foliofont\thepage}\hfill{\lrhfont\rightmark}\hfill {\@volname \@jvol}}%
      \def\@oddhead{\@articlebar}%
      %\def\@oddhead{\GRIDRULES\llap{\hrule height .5pt width \textwidth}{\lrhfont\thepage}\hfill{\lrhfont\rightmark}\hfill {\lrhfont \@volname \@jvol}}%
      \def\@oddfoot{}
      \def\@evenfoot{}
      \let\@mkboth\markboth
      \def\sectionmark##1{\markboth{##1}{##1}}%
       }
\else
  \def\ps@headings{%
    \let\@oddfoot\@empty
    \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
    \let\@mkboth\markboth
    \def\sectionmark##1{\markright {##1}}}
\fi
\def\crpageinfo#1{\gdef\@crpageinfo{#1}}
\crpageinfo{\hfil{\foliofont\thepage}}
\def\crpageinfoverso#1{\gdef\@crpageinfoverso{#1}}
\crpageinfoverso{{\foliofont\thepage}\hfil\copyrightfont
\copyright\ The Author \@copyrightyear. Published by Oxford
University Press. All rights reserved. For Permissions, please
email: journals.permissions@oup.com}
%\def\ps@plain{\let\@mkboth\@gobbletwo
%     \def\@oddhead{\GRIDRULES%
%     \raisebox{-13pt}{\smash{\hbox to \textwidth{\hskip\jnameskip{\jnamefont{\@jname}}\hfill\hfill{}}}}
%     %
%     \raisebox{-28pt}{\llap{\hskip-35pc\jnamefont Copyright\copyright\ Socity of System Biologyst\hskip35pc\hfill}}}%
%     \def\@evenhead{\GRIDRULES%
%     \raisebox{-13pt}{\smash{\hbox to \textwidth{\hskip\jnameskip{\jnamefont{\@jname}}\hfill{\ifx\@articletype\articletypecheck\articletypesmallfont\else\articletypefont\fi\@articletype}\hfill{\catchlineinfofont\catchlineinfo}}}}%
%     \llap{\rule[-26pt]{\textwidth}{1pt}}}%
%      \def\@oddfoot{\rlap{\rule[12pt]{\textwidth}{1pt}}\copyrightfont \@crpageinfo}
%      \def\@evenfoot{\rlap{\rule[12pt]{\textwidth}{1pt}}\@crpageinfoverso}
%%      \def\@evenfoot{\rlap{\rule[12pt]{\textwidth}{1pt}}{\foliofont\thepage}\hfil\copyrightfont \copyright\ The Author 2016. Published by Oxford University Press. All rights reserved. For Permissions, please email: journals.permissions@oxfordjournals.org}
%     }

\newdimen\cropenskip\cropenskip=15pt

\def\cropeninfo#1{\gdef\@cropeninfo{#1}}

\def\jshort#1{\gdef\@jshort{#1}}
\jshort{sys}

\def\pvtitle#1{\gdef\@pvtitle{#1}}
\pvtitle{Point of View}


\newdimen\PVskip\PVskip=-55.2pt


\if@mstype
\cropeninfo{{\fontsize{6}{7}\selectfont\hspace*{-7.3pt}
\begin{tabular}[t]{l}
\\[\PVskip]
\begin{boxedminipage}{\textwidth}{\vspace*{14.4pt}\begin{center}
\titlefont \@pvtitle
\end{center}\vspace*{4pt}}
\end{boxedminipage}\\[34pt]
{\it Syst. Biol.} \@jvolume(\@jissue):\@firstpage\@pageinfo, \@pubyear\\
\copyright\ The Author(s) 2016. Published by Oxford University Press, on behalf of the Society of Systematic Biologists. All rights reserved.\\
For Permissions, please email: journals.permissions@oup.com\\
DOI:10.1093/sysbio/\@jshort\@artid\\
\accessfont{\@access}
\end{tabular}}}
\else
\cropeninfo{{\fontsize{7}{9}\selectfont\hspace*{-6.3pt}\begin{tabular}[t]{l}
\\[-30.4pt]
\if@bookreview\\[48.2pt]\else\\[-15.4pt]\fi
%{\it Syst. Biol.gdsgfs} \@jvolume(\@jissue):\@firstpage\@pageinfo, \@pubyear\\
\copyright\ The Author 2016. Published by Oxford University Press on behalf of the Society for Molecular Biology and Evolution. All rights reserved.\\ For permissions, please email: journals.permissions@oup.com\\
\rule[4pt]{\textwidth}{.5pt}\\
{\fontsize{10}{10}\selectfont{\it Mol. Biol. Evol.} \@jvolume(\@jissue):\@firstpage\@pageinfo\quad doi:10.1093/molbev/\@jshort\@artid}\quad \accessfont{\@access\hfill{\foliofont \thepage}}\\[38pt]
\end{tabular}}}
\fi




\def\ps@plain{\let\@mkboth\@gobbletwo
    % \def\@oddhead{\GRIDRULES\@cropeninfo}%
    % \def\@evenhead{\GRIDRULES\@cropeninfo}%
     \def\@oddhead{\GRIDRULES\@articlebar}%
     \def\@evenhead{\GRIDRULES}%
      \def\@oddfoot{\@cropeninfo\par}
      \def\@evenfoot{\@cropeninfo\hfill\thepage}
     }

%\def\crpageinfo#1{\gdef\@crpageinfo{#1}}
%\crpageinfo{\copyright\ The Author \@copyrightyear. Published by
%Oxford University Press. All rights reserved. For Permissions,
%please email:
%journals.permissions@oxfordjournals.org\hfil{\foliofont\thepage}}
%\def\ps@plain{\let\@mkboth\@gobbletwo
%     \def\@oddhead{\GRIDRULES%
%     \raisebox{-13pt}{\smash{\hbox to \textwidth{\hskip\jnameskip{\jnamefont{\@jname}}\hfill{\ifx\@articletype\articletypecheck\articletypesmallfont\else\articletypefont\fi\@articletype}\hfill{\catchlineinfofont\catchlineinfo}}}}%
%     \llap{\rule[-26pt]{\textwidth}{1pt}}}%
%     \def\@evenhead{\GRIDRULES%
%     \raisebox{-13pt}{\smash{\hbox to \textwidth{\hskip\jnameskip{\jnamefont{\@jname}}\hfill{\ifx\@articletype\articletypecheck\articletypesmallfont\else\articletypefont\fi\@articletype}\hfill{\catchlineinfofont\catchlineinfo}}}}%
%     \llap{\rule[-26pt]{\textwidth}{1pt}}}%
%      \def\@oddfoot{\rlap{\rule[12pt]{\textwidth}{1pt}}\copyrightfont \@crpageinfo}
%      \def\@evenfoot{\rlap{\rule[12pt]{\textwidth}{1pt}}{\foliofont\thepage}\hfil\copyrightfont \copyright\ The Author 2008. Published by Oxford University Press. All rights reserved. For Permissions, please email: journals.permissions@oxfordjournals.org}
%     }
%
\def\ps@front{\let\@mkboth\@gobbletwo
     \def\@oddhead{\GRIDRULES%
     \raisebox{-9pt}{\smash{\hbox to \textwidth{\hskip\jnameskip{\jnamefont{\@jname}}\hfill{\ifx\@articletype\articletypecheck\articletypesmallfont\else\articletypefont\fi\@articletype}\hfill{\catchlineinfofont\catchlineinfo}}}}%
     \llap{\rule[-16pt]{\textwidth}{1pt}}}%
     \def\@evenhead{\GRIDRULES%
     \raisebox{-9pt}{\smash{\hbox to \textwidth{\hskip\jnameskip{\jnamefont{\@jname}}\hfill{\ifx\@articletype\articletypecheck\articletypesmallfont\else\articletypefont\fi\@articletype}\hfill{\catchlineinfofont\catchlineinfo}}}}%
     \llap{\rule[-16pt]{\textwidth}{1pt}}}%
      \def\@oddfoot{\rlap{\rule[12pt]{\textwidth}{1pt}}\copyrightfont }
      \def\@evenfoot{\rlap{\rule[12pt]{\textwidth}{1pt}}\copyrightfont }
     }
\def\ps@contents{%
      %\let\@oddfoot\@empty
      %\let\@evenfoot\@empty
      \def\@evenhead{\GRIDRULES\rlap{\rule[-10pt]{\textwidth}{.5pt}}{\lrhfont\leftmark}\hfil}%
      \def\@oddhead{\GRIDRULES\rlap{\rule[-10pt]{\textwidth}{.5pt}}\hfil{\rrhfont\rightmark}}%
      \def\@oddfoot{\rlap{\rule[12pt]{\textwidth}{1pt}}}
      \def\@evenfoot{\rlap{\rule[12pt]{\textwidth}{1pt}}}
    \let\@mkboth\@gobbletwo
    \let\chaptermark\@gobble
    \let\sectionmark\@gobble
    \def\titlemark##1{\markboth{##1}{##1}}%
    \def\authormark##1{\gdef\leftmark{##1}}%
}

\def\jtitlemark#1{\gdef\@jtitlemark{#1}}

\jtitlemark{SYSTEMATIC BIOLOGY}

\font\rrrhfont=cmr10 at 18pt

\def\articlebar#1{\gdef\@articlebar{#1}}

\articlebar{\begin{picture}(0,0)(0,0)
\rotatebox{-90}{\put(415,514.5){\rule{88pt}{43pt}}
\put(415,525){\hbox to 88pt{\hfill\fontsize{20}{20}\bfseries\selectfont\textcolor{white}{\@mstype}\hfill}}}
\end{picture}}

\def\aarticlebar#1{\gdef\@aarticlebar{#1}}

\aarticlebar{\begin{picture}(0,0)(0,0)
\rotatebox{-90}{\put(520,500){\color{grayfifty}\rule{88pt}{43pt}}
\put(520,510){\hbox to 88pt{\hfill\fontsize{20}{20}\bfseries\selectfont\textcolor{white}{Fast Track}\hfill}}}
\end{picture}}

\def\ps@myheadings{%
      %\let\@oddfoot\@empty
      %\let\@evenfoot\@empty
      \def\@evenhead{\GRIDRULES\rlap{\rule[-6pt]{\textwidth}{.5pt}}{\lrhfont\leftmark\enskip$\cdot$\enskip \lrhifont{doi:10.1093/molbev/\@jshort\@artid}}\hfil{}\hfil{\rrrhfont {M\hskip0.8ptB\hskip0.3ptE}}}%
      \def\@oddhead{\GRIDRULES\rlap{\rule[-6pt]{\textwidth}{.5pt}}{\rrhfont\rightmark\enskip$\cdot$\enskip \lrhifont{doi:10.1093/molbev/\@jshort\@artid}}\hfil{}\hfil{\rrrhfont {M\hskip0.8ptB\hskip0.3ptE}}}%
      \def\@oddfoot{\hfill\foliofont\thepage}
      \def\@evenfoot{\foliofont\thepage}
    \let\@mkboth\@gobbletwo
    \let\chaptermark\@gobble
    \let\sectionmark\@gobble
    \def\titlemark##1{\markboth{##1}{##1}}%
    \def\authormark##1{\gdef\leftmark{##1}}%
}
%
\def\ps@empty{\let\@mkboth\@gobbletwo%
  \def\@oddhead{\GRIDRULES}%
  \let\@oddfoot\@empty%
  \def\@evenhead{\GRIDRULES}%
  \let\@evenfoot\@empty}
%
%%%%%%%%%%%%%%%%%%%% ** COLOR DEFINITIONS ** %%%%%%%%%%%%%%%%%%%
%
%\usepackage[dviwindo,pantonenames]{colorx}
%\definecolor{PMS300C}{named}{Pantone 300 C}
%03-12-13%\definecolor{C94M89Y0K0}{named}{Pantone C94M89Y0K0}
%03-12-13%\definecolor{C60M80Y0K50}{named}{Pantone C60M80Y0K50}
%B%\definecolor{PMS300C40}{named}{Pantone 300 40C}
%B%\definecolor{myblue}{named}{Blue}
%B%\definecolor{snamecolor}{cmyk}{0,1,0,0}
%B%\definecolor{tint}{cmyk}{0,0,0,.1}
%B%\definecolor{blue}{rgb}{0,0,1}
%B%\definecolor{lightblue}{rgb}{0,0,.6}
\definecolor{grayten}{cmyk}{0,0,0,.1}
%B%\definecolor{graytwenty}{cmyk}{0,0,0,.2}
%B%\definecolor{grayfifteen}{cmyk}{0,0,0,.15}
%B%\definecolor{graytwentyfive}{cmyk}{0,0,0,.25}
%B%\definecolor{graythirty}{cmyk}{0,0,0,.3}
\definecolor{grayfifty}{cmyk}{0,0,0,.5}
%B%%
%%%%%%%%%%%%%%%%%%%% ** LOCAL FONT DEFINITIONS ** %%%%%%%%%%%%%%%%%%%
%
% Running head, foot, pageno & other info.
\def\lrhfont{\fontfamily{\rmdefault}\fontseries{m}\fontshape{n}\fontsize{10}{10}\selectfont}
\def\rheadfont{\fontfamily{\rmdefault}\fontseries{m}\fontshape{n}\fontsize{18}{18}\selectfont}
\def\lrhifont{\fontfamily{\rmdefault}\fontseries{m}\fontshape{n}\fontsize{10}{10}\selectfont}
\def\rrhfont{\fontfamily{\rmdefault}\fontseries{m}\fontshape{n}\fontsize{10}{10}\selectfont}
\def\rrffont{\fontfamily{\rmdefault}\fontseries{b}\fontshape{it}\fontsize{7}{10}\selectfont}
\def\lrffont{\fontfamily{\rmdefault}\fontseries{b}\fontshape{it}\fontsize{7}{10}\selectfont}
\def\foliofont{\reset@font\fontfamily{\rmdefault}\fontseries{m}\fontshape{n}\fontsize{10}{12}\selectfont}
\def\opheadfootfont{\fontfamily{\rmdefault}\fontsize{8}{8}\selectfont}
% Header
\def\chapwordfont{\fontfamily{\sfdefault}\fontsize{14}{12}\fontseries{x}\selectfont\leftskip0pt\rightskip0pt plus1fill}
\def\chapnofont{\fontfamily{\sfdefault}\fontsize{42}{42}\fontseries{cn}\fontshape{n}\selectfont}
\def\chapterfont{\fontfamily{\sfdefault}\fontsize{51.5}{54}\fontseries{cn}\fontshape{n}\selectfont\leftskip0pt plus 1fill\rightskip0pc}
\def\subchapterfont{\fontfamily{\rmdefault}\fontsize{14}{17}\fontseries{ds}\fontshape{it}\selectfont\leftskip0pt plus 1fill\rightskip0pt plus 1fill}
\def\unnumchapterfont{\fontfamily{\rmdefault}\fontsize{12}{12}\fontseries{b}\fontshape{n}\selectfont\leftskip0pt\rightskip0pt plus 1fill}
\def\titlefont{\fontfamily{\rmdefault}\fontsize{20}{22}\fontseries{b}\selectfont\leftskip0pt\rightskip0pt plus1fill}
\def\subtitlefont{\fontfamily{\rmdefault}\fontsize{14}{18}\fontseries{m}\selectfont\leftskip0pt\rightskip0pt plus1fill}
\def\authorfont{\reset@font\fontfamily{\rmdefault}\fontsize{12.5}{14}\fontseries{m}\selectfont\leftskip0pt\rightskip0pt plus1fill}
\def\fmauthorfont{\fontfamily{\sfdefault}\fontsize{12}{18}\fontseries{m}\fontshape{n}\selectfont\leftskip36pt\rightskip0pt plus1fill}%\reset@font
\def\addressfont{\reset@font\fontfamily{\rmdefault}\fontshape{m}\fontsize{10}{14}\selectfont\leftskip0pt\rightskip0pt plus1fill}
\def\abstractfont{\fontfamily{\rmdefault}\fontsize{10.5}{20}\selectfont\leftskip0pt\rightskip0pt}
\def\abstractheadfont{\fontfamily{\rmdefault}\fontseries{b}\fontshape{n}\fontsize{13}{14}\selectfont}
\def\keywordheadfont{\fontfamily{\sfdefault}\fontseries{b}\fontsize{10.5}{12}\selectfont}
\def\keywordfont{\fontfamily{\rmdefault}\fontsize{10.5}{12}\selectfont\leftskip0pt\rightskip0pt}
\def\historyfont{\fontfamily{\rmdefault}\fontsize{8}{10}\fontshape{it}\selectfont\leftskip0pt\rightskip0pt plus1fill\centering}
%\def\accessfont{\fontfamily{\rmdefault}\fontsize{6}{8}\fontseries{l}\selectfont\leftskip0pt\rightskip0pt plus1fill}
\def\accessfont{\fontfamily{\rmdefault}\fontsize{10}{10}\selectfont\fontshape{rmdefault}\selectfont}
\def\corespfont{\fontfamily{\rmdefault}\fontsize{10}{14}\fontshape{m}\selectfont\leftskip0pt\rightskip0pt plus1fill}
\def\editorfont{\fontfamily{\rmdefault}\fontsize{10}{14}\fontshape{m}\selectfont\leftskip0pt\rightskip0pt plus1fill}
% Sections
\def\sectionfont{\fontfamily{\rmdefault}\fontsize{13}{15}\fontseries{b}\selectfont\leftskip0pt\rightskip0pt plus1fill}
\def\subsectionfont{\fontfamily{\rmdefault}\fontsize{11.5}{13}\fontshape{n}\selectfont\raggedright}
\def\subsubsectionfont{\fontfamily{\rmdefault}\fontshape{it}\fontsize{10.5}{12}\selectfont\raggedright}
\def\subsubparagraphfont{\fontfamily{\rmdefault}\fontshape{rm}\selectfont\raggedright}
\def\paragraphfont{\fontfamily{\rmdefault}\fontshape{it}\selectfont\raggedright}
\def\subparagraphfont{\fontfamily{\rmdefault}\fontsize{8}{10}\fontseries{b}\selectfont}
\def\listheadfont{\normalfont\normalsize}
% Floats
\def\timesc{\fontsize{7.5}{12}\selectfont}
\def\timescm{\fontsize{6}{10}\selectfont}
\def\figurenumberfont{\reset@font\fontsize{9}{9}\fontseries{b}\selectfont}%
\def\figurecaptionfont{\fontsize{9}{9}\selectfont}
%
\def\tablenumberfont{\reset@font\fontfamily{\rmdefault}\fontsize{9}{9}\fontseries{b}\selectfont}
\def\tablecaptionfont{\fontfamily{\rmdefault}\fontsize{9}{10}\selectfont\leftskip\tablecapleftskip\rightskip\tablecaprightskip}
\def\tablebodyfont{\fontfamily{\sfdefault}\fontseries{b}\fontsize{8.5}{12}\selectfont\leftskip\tableleftskip\rightskip\tablerightskip}
\def\tablefootfont{\fontfamily{\rmdefault}\fontsize{7.5}{9}\selectfont\leftskip\tablefootleftskip \rightskip\tablefootrightskip}
\def\TCHfont{\fontfamily{\rmdefault}\fontsize{8}{9}\selectfont}
\def\SCHfont{\fontfamily{\rmdefault}\fontsize{8}{9}\selectfont}
% List
\def\listfont{\normalfont}
\def\listlabelfont{\normalfont}
\def\mclistfont{\fontsize{10}{13}\selectfont}
% Biliography
\def\bibfont{\fontsize{8}{18}\selectfont}%9
% Index
\def\indexfont{\fontsize{9}{11}\selectfont}
% OTHERS
\def\extractfont{}%\fontsize{10}{12}\selectfont
\def\thmfont{\fontfamily{\rmdefault}\fontsize{10}{20}\fontshape{it}\selectfont}
\def\quotefont{\fontsize{10}{11}\selectfont}
\def\symbolsfont{\fontsize{10}{13}\selectfont}
\def\ackfont{\fontsize{10}{12}\selectfont}

%
%%%%%%%%%%%%%%%%%%%% ** JOURNAL HEADER COMMANDS ** %%%%%%%%%%%%%%%%%%%
%
%\def\jyear#1{\gdef\@jyear{20#1}\gdef\@shortyear{#1}}

\newskip\jnameskip
\def\jnamefont{\fontfamily{\rmdefault}\fontsize{6}{8}\selectfont\fontshape{it}\selectfont}
\def\jname#1{\gdef\@jname{#1}}\jname{Syst. Biol.}
%
\def\articletypefont{\fontfamily{\rmdefault}\fontsize{19}{21}\selectfont\fontseries{b}\fontshape{it}\selectfont}
\def\articletypesmallfont{\fontfamily{\rmdefault}\fontsize{17}{21}\selectfont\fontseries{b}\fontshape{it}\selectfont}
\def\articletype#1{\gdef\@articletype{#1}}
\articletype{ORIGINAL PAPER}
\def\articletypecheck{APPLICATIONS NOTE}
%
\def\articlesubsubtypefont{\fontfamily{\rmdefault}\fontsize{21}{23}\selectfont\fontshape{it}\selectfont\bfseries}
\def\articlesubsubtype#1{\gdef\@articlesubsubtype{#1}}
\articlesubsubtype{}
%
\def\articlesubtypefont{\fontfamily{\rmdefault}\fontsize{12}{14}\selectfont\fontshape{it}\selectfont}
\def\articlesubtype#1{\gdef\@articlesubtype{#1}}
\articlesubtype{Sequence analysis}
%
\def\jvolume#1{\gdef\@jvolume{#1}}%\jvolume{00}
\def\jvol#1{\gdef\@jvol{#1}}%\jvol{00}
\def\jissue#1{\gdef\@jissue{#1}}%\jissue{00}
\def\pubyear#1{\def\@pubyear{#1}}\pubyear{0000}

\def\copyrightfont{\fontfamily{\rmdefault}\fontsize{7}{9}\selectfont\fontshape{m}\selectfont}
\def\copyrightyear#1{\def\@copyrightyear{#1}}
\def\artid#1{\gdef\@artid{#1}}\artid{000}
\def\shortjname#1{\gdef\@shortjname{#1}}\shortjname{btm}
%
\def\pages#1{\gdef\@pages{#1}}
\pages{pages}
\def\pageinfo#1{\gdef\@pageinfo{#1}}
\pageinfo{{--}\@lastpage}
\def\catchlineinfofont{\fontfamily{\rmdefault}\fontsize{8}{10}\selectfont\fontshape{it}\selectfont}
\def\catchlineinfo{\begin{tabular}[b]{@{}r@{}}
\ifx\@jvol\undefined\else Vol.\ \@jvol\ \fi\ifx\@jissue\undefined\else no.\ \@jissue\ \fi\@pubyear,\ \@pages\ \@firstpage\@pageinfo\\%\ifnum\@lastpage=\@firstpage\relax\else--\@lastpage\fi
%\fi\@pubyear,\ pages\ \@firstpage--\@lastpage\\
doi:10.1093/sysbio/\@shortjname\@artid
\end{tabular}}
%
\long\def\title{\@ifnextchar[{\short@title}{\@@title}}
\def\short@title[#1]{\titlemark{#1}\@@@title}
\def\@@title#1{\authormark{#1}\@@@title{#1}}
\long\def\@@@title#1{\gdef\@title{#1}\addcontentsline{toc}{article}{#1}}
%
\def\subtitle#1{\gdef\@subtitle{{#1}}}
%
\long\def\author{\@ifnextchar[{\short@uthor}{\@uthor}}
\def\short@uthor[#1]{\authormark{#1}\@@author}
\def\@uthor#1{\authormark{#1}\@@author{#1}}
\long\def\@@author#1{\gdef\@author{#1}}
\long\def\surname#1{\if@surnamecolor{\textcolor{snamecolor}{#1}}\else{{#1}}\fi}
\def\cename#1{\gdef\@cename{#1}}
\cename{\relax}
\def\copyeditor{\if@CE{\protect{\slugfont Copyedited by: \@cename}}\else\relax\fi}
\def\mstype#1{\gdef\@mstype{#1}}
\mstype{\relax}
\def\msc{\if@MSC{\protect{\slugfont MANUSCRIPT CATEGORY:\quad \@mstype}}\else\relax\fi}
%
\def\address#1{\gdef\@address{#1}}
%
\def\abstract#1{\gdef\@abstract{\parindent=10pt #1\par\nobreak}}
%
\def\history#1{\gdef\@history{#1}}\history{}
\def\access#1{\gdef\@access{#1}}\access{}
\def\coresp#1{\gdef\@coresp{$^{\boldsymbol\ast}$\textbf{Corresponding author:\ }#1}}\coresp{}
\def\editor#1{\gdef\@editor{\textbf{Associate Editor:\ }#1}}\editor{}
\def\datade#1{\gdef\@datade{#1}}\datade{}

%
\def\rec#1{\unskip Received #1}%
\def\acc#1{\unskip Accepted #1}%

\def\keyword#1{\def\@key{{\keywordheadfont \keywordname}:\keywordfont~#1}}
\gdef\email#1{E-mail: #1}

%Abstract macro been modified by Ganga on 23rd June 2010
%\newenvironment{abstract}[1][\abstractname]{\begingroup%
%  \abstractfont%
%  \let\section\subparagraph%
%  \def\@periodaftersection{}%
%  \ifx#1\empty{\fontseries{b}\selectfont \abstractname}\else{\fontseries{b}\selectfont #1}\fi\@abstract}
%{\par\vskip6pt\noindent\@history\endgroup}
%\bigskip\@afterheading\@afterindenttrue

%
\def\mystar{\ifx\@webpdf\undefined
            $\ast$
            \else%
            \hyperlink{\csname fn\@roman\c@footnote\endcsname}{\textsuperscript{$\ast$}}%
            \fi}
          %
\def\@fnsymbol#1{\ensuremath{\ifcase#1\or {\ast} \or {\dagger}\or \ddagger\or
   \mathsection\or \mathparagraph\or \|\or **\or \dagger\dagger
   \or \ddagger\ddagger \else\@ctrerr\fi}}
%
\def\@myfnsymbol#1{\ensuremath{\ifcase#1\or ^{\ast} \or ^{\dagger}\or ^{\ddagger}\or
   \mathsection\or \mathparagraph\or \|\or **\or \dagger\dagger
   \or \ddagger\ddagger \else\@ctrerr\fi}}
%
%%%%%%%%%%%%%%%%%%%% ** ARTICLE HEADER ** %%%%%%%%%%%%%%%%%%%
%
\newcommand\maketitle{\par
  \begingroup
  \renewcommand\thefootnote{\@myfnsymbol\c@footnote}%
    \long\def\@makefntext##1{\footnotesize%
    \if\c@footnote=1%
       \ifx\@webpdf\undefined%
          {\mystar}%
       \else%
          \hypertarget{\csname fn\@roman\c@footnote\endcsname}{\mystar}%
       \fi%
    \else%
       \ifx\@webpdf\undefined%
          {{\@thefnmark}}%
       \else%
          \hypertarget{\csname fn\@roman\c@footnote\endcsname}{{\textsuperscript{\@thefnmark}}}%
       \fi##1}%
    \def\@makefnmark{%
    \if\c@footnote=1%
       \ifx\@webpdf\undefined%
          {\mystar}%
       \else%
          \hyperlink{\csname fn\@roman\c@footnote\endcsname}{\mystar}%
       \fi%
    \else%
       \ifx\@webpdf\undefined%
          {{\@thefnmark}}%
       \else%
          \hyperlink{\csname fn\@roman\c@footnote\endcsname}{{\@thefnmark}}%
       \fi%
    \fi}%
    \if@twocolumn%
      \ifnum \col@number=\@ne%
        \@maketitle%
      \else%
        \twocolumn[\@maketitle]%
      \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \thispagestyle{plain}\@thanks
      \endgroup
  \@afterindenttrue
  \setcounter{footnote}{0}%
  \global\let\@chapter\@empty\global\let\type\relax
  \global\let\@title\@empty\global\let\title\relax
  \global\let\@subtitle\@empty\global\let\subtitle\relax %%subtitle
  \global\let\@author\@empty\global\let\author\relax
  \global\let\@thanks\@empty\global\let\thanks\relax
  \global\let\@address\@empty\global\let\address\relax
  \global\let\affD\relax
  \global\let\@Date\@empty\global\let\Date\relax
  \global\let\@abstract\@empty\global\let\abstract\relax
  \global\let\@key\@empty\global\let\key\relax
  \global\let\@tableofcontents\@empty\global\let\tableofcontents\relax
  \global\let\maketitle\relax\global\let\@maketitle\relax
\@afterheading
\enlargethispage{-36pt}}
%
\newdimen\titleskip\if@mstype\titleskip=26pt\else\titleskip=-.8pt\fi
\newdimen\authorskip\if@mstype\authorskip=-20pt\else\authorskip=-12.5pt\fi
\newdimen\ruleskip\ruleskip=4pt
\newdimen\ruleht\ruleht=1pt
\newskip\belowmktitskip\belowmktitskip=3.05pt
\newskip\alignmkadj
\def\@maketitle{%
 \begin{flushleft}%
  \let \footnote\thanks
 \ifx\@articlesubtype\empty\else
    \vspace*{-32pt}%

 \fi%
 \ifx\@title\empty\else
    \vspace*{\titleskip}%
    {\titlefont\@title\par}
 \fi%
 \ifx\@subtitle\undefined\else
     \vskip5.7pt
    {\subtitlefont\@subtitle\par}
 \fi%
 \ifx\@author\empty\else
    \vskip\authorskip%
    {\authorfont\@author\par}
 \fi%
 \ifx\@address\empty\else
    \vspace*{1pt}%
        {\addressfont\@address\par}
        %sen% \fi%
 \ifx\@coresp\empty\else
    \vspace*{0pt}%
    {\corespfont\@coresp\par}
 \fi%
%B%\ifx\@history\empty\else
%B% \vspace*{7pt}%
%B%{\historyfont\@history\par}
%B%\fi%
 \ifx\@editor\empty\else
    \vspace*{-1pt}%
    {\editorfont\@editor\par}
 \fi%
\ifx\@datade\empty\else
    \vspace*{0pt}%
    {\corespfont\@datade\par}
 \fi%

%sen% \ifx\@access\empty
%sen%    \vspace*{8.5pt}%
%sen%    \else
%sen%    \vspace*{3pt}%
%sen%    {\accessfont\@access\par}

\fi%
\if@mstype\else
 \ifx\@abstract\empty\else
    \vspace*{14.8pt}%
    {\abstractname\vskip6pt\par\noindent\abstractfont\@abstract\par\vspace*{22pt}}
\fi\fi
 \ifx\@key\empty\else
    \vspace*{-9pt}%
    {\editorfont\@key\par}
    \vspace*{15pt}
 \fi%
%\vskip-13.2pt
%\vskip-\ruleskip\vrule height\ruleht depth0pt width\textwidth\par
  \advance\belowmktitskip\alignmkadj
  \vskip\belowmktitskip
\end{flushleft}%
\par
}
%
%
\newcommand\frontmatter{%
  %\cleardoublepage
  %\@mainmatterfalse
  \twocolumn
  \pagenumbering{roman}
  \pagestyle{fmbm}}
\newcommand\mainmatter{%
    \cleardoublepage
  \@mainmattertrue
  \pagenumbering{arabic}}
\newcommand\backmatter{%
  \if@openright
    \cleardoublepage
  \else
    \clearpage
  \fi
  \@mainmatterfalse}
\newcommand\chapter{\if@openright\cleardoublepage\else\clearpage\fi
                    \thispagestyle{plain}%
                    \global\@topnum\z@
                    \@afterindentfalse
                    \secdef\@chapter\@schapter}
\def\@chapter[#1]#2{\let\makemytitle\maketitle
                      \ifnum \c@secnumdepth >\m@ne
                       \if@mainmatter
                         \refstepcounter{chapter}%
                         \typeout{\@chapapp\space\thechapter.}%
                         \addcontentsline{toc}{chapter}%
                                   {\protect\numberline{\thechapter}#1}%
                       \else
                         \addcontentsline{toc}{chapter}{#1}%
                       \fi
                    \else
                      \addcontentsline{toc}{chapter}{#1}%
                    \fi
                    %\bookmark{Bounded Rationality and Industrial Organization}
                    %\bookmark{\ifnum\c@chapter > 11 Discussion \else \ifnum\c@chapter > 8 Reference Dependence\else  \ifnum\c@chapter > 5 Responding to Market Complexity \else  \ifnum\c@chapter > 1 Anticipating Future Preferences\else  Introduction\fi  \fi \fi\fi}
                    \chaptermark{#1}
                    \addtocontents{lof}{\protect\addvspace{10\p@}}%
                    \addtocontents{lot}{\protect\addvspace{10\p@}}%
                    \gdef\@chaptertitle{{{#2}}}}
%
%
\newcommand{\makemytitle}{
                    \if@twocolumn
                      \@topnewpage[\@makechapterhead{\@chaptertitle}]%
                    \else
                      \@makechapterhead{\@chaptertitle}%
                    \fi%
                    \@afterheading%
                    \gdef\@author{}%
                    }
%
%\let\makemytitle\maketitle
%
\newskip\alignadj
\newskip\aboveheaderskip\aboveheaderskip=0pt
\newskip\belowheaderskip\belowheaderskip=35pt
\newskip\abovefmheaderskip\abovefmheaderskip=0pt
\newskip\belowfmheaderskip\belowfmheaderskip=59.3pt

\newskip\abovecnskip\abovecnskip=24\p@
\newskip\belowcnskip\belowcnskip=0\p@
\newskip\abovectskip\abovectskip=29pt
\newskip\belowctskip\belowctskip=-3\p@

\newskip\aboveauskip\aboveauskip=0pt%
\newskip\belowauskip\belowauskip=0\p@

\newskip\aboveaffskip\aboveaffskip=0\p@
\newskip\belowaffskip\belowaffskip=0\p@
%
\newskip\abovefmcnskip\abovefmcnskip=2\p@
\newskip\belowfmcnskip\belowfmcnskip=0\p@
\newskip\abovefmctskip\abovefmctskip=-26\p@
\newskip\belowfmctskip\belowfmctskip=10.3\p@
%
\newif\ifminspace\global\minspacefalse
\def\minimumspace{\global\minspacetrue}
%
\newskip\firstlinefixedpos\firstlinefixedpos=0pt
%
\def\myhrulefill{\leavevmode\leaders\hrule \@height 1pt\hfill\kern\z@}
%
\def\subtitle#1{\gdef\@subtitle{#1\vphantom{ph}}}\subtitle{}%
\def\source#1{\gdef\@source{\noindent#1}}
%
\def\chaptoc#1{\gdef\@chaptoc{#1}}\chaptoc{}%
%
\newbox\chapnobox
\newdimen\chapwdwidth%
\def\rulewidth#1{\def\@rulewidth{#1}}
\newdimen\chapnowd
%
\def\@makechapterhead#1{%
\setbox\chapnobox=\hbox{\chapnofont{\@chapapp\hskip6.4mm\Thechapter}}
\chapnowd\wd\chapnobox
  \ifdim\firstlinefixedpos>0pt \vbox to \firstlinefixedpos\else \vbox\fi{\kern\aboveheaderskip%
    \ifnum \c@secnumdepth >\m@ne
      \if@mainmatter
       \vskip\abovecnskip
       \noindent{\fboxsep=10pt\colorbox{cyan}{\hbox to 184mm{\hfill{\chapnofont\textcolor{white}{\@chapapp\space\thechapter}\hskip25mm}}\par}}%
        %\vskip\belowcnskip
       \interlinepenalty\@M
       \vskip\abovectskip
       {{\noindent\chapterfont{#1}\vphantom{ph}\par}}
       \vskip0pt
       \ifx\@subtitle\@undefined
       %\vskip-5pt
       \belowctskip=0pt
       \else
%      \vskip19.7pt
       \vskip0pt
       \noindent{\subchapterfont\begin{minipage}[t]{\textwidth}\centering\@subtitle\par\end{minipage}}
       %\vskip17pt
       \fi
       \vskip\belowctskip
       \ifx\@author\@undefined
       \else
       \vskip\belowctskip
       \vskip0pt
%      \vskip31.6pt
       \noindent{\authorfont\begin{minipage}[t]{\textwidth}\centering\@author\par\end{minipage}}
       \fi
  \advance\belowheaderskip\alignadj
  \vskip\belowheaderskip
  \fi\fi
  }\vskip \z@ plus 2pt}
%
\def\@schapter#1{\gdef\@chaptertitle{#1}\let\makemytitle\makestitle\chaptermark{#1}}

\newcommand{\makestitle}{%
                 \if@twocolumn%
                   \@topnewpage[\@makeschapterhead{\@chaptertitle}]%
                 \else%
                   \@makeschapterhead{\@chaptertitle}%
                 \fi%
                 \@afterheading%
                 \let\@author\relax%
                 }
%
\newdimen\chrulesp
\chrulesp=-4pt
\def\@makeschapterhead#1{\thispagestyle{plain}%\firstlinefixedpos=214.5pt%
  \ifdim\firstlinefixedpos>0pt \vbox to \firstlinefixedpos\else \vbox\fi{\kern\aboveheaderskip%
    \ifnum \c@secnumdepth >\m@ne
      \if@mainmatter
        \vskip\abovefmcnskip
        \vskip-1pt
        %{\noindent\raise 6pt\hbox to 6pc{\chapwordfont{\vphantom{\@chapapp}}}\hskip1pc\hbox to \textwidth{\chapnofont\noindent\vphantom{\thechapter}}}%
        %\vskip\belowcnskip
        %\noindent\textcolor{graytwenty}{\rule{3pc}{2pt}}%
      \fi
    \fi
    \interlinepenalty\@M
      \vskip\abovefmctskip
      \noindent{\fboxrule=.5pt\fboxsep=18pt\fbox{\unnumchapterfont\begin{minipage}[t]{465pt}\centering{{{#1}}}\vspace{-3pt}\vphantom{ph}\par\end{minipage}}}
      \vskip\belowfmctskip
      \ifx\@author\empty\else
        \vskip\aboveauskip
        {\noindent{\fmauthorfont\@author\par}}
        \vskip\belowauskip
      \fi
  \advance\belowfmheaderskip\alignadj
  \vskip\belowfmheaderskip
  }\vskip 0\p@}
%
\newcommand{\myarticletype}[1]{{\hskip10pt\normalsize\bfseries\selectfont{#1}\vphantom{phy}}\par\nobreak}
\newcommand{\myarticleau}[1]{{\rule[2pt]{10pt}{.75pt}\normalsize\normalfont\selectfont{#1}\vphantom{phy}}}
\newcommand{\myarticleaff}[1]{{\normalsize\normalfont\selectfont{#1}\vphantom{phy}}\par\nobreak\vskip.2\baselineskip}
\newcommand{\myhistory}[1]{\par\vskip2\baselineskip\noindent{\footnotesize\normalfont\itshape\selectfont{#1}\vphantom{phy}}\par\nobreak}
\newcommand{\myconrrigendum}[1]{\par\vskip2\baselineskip\leftskip10pt\noindent{\fontsize{10}{12}\selectfont{#1}\vphantom{phy}}\par\nobreak}
\newcommand{\myconrrigendumtext}[1]{\par\vskip2\baselineskip\leftskip0pt{\fontsize{10}{12}\selectfont{#1}\vphantom{phy}}\par\nobreak}

%
%%%%%%%%%%%%%%%%%%%% ** DROPCAP ** %%%%%%%%%%%%%%%%%%%
%
\def\dropfont{\fontsize{26}{26}\selectfont\bfseries}
\newcount\droplines\droplines=2
%
\newdimen\dropsep\dropsep=2pt
\newskip\capvskip\capvskip=5pt
\newskip\caphskip\caphskip=0pt
%
\newdimen\capwd
\newdimen\droplinewd
%
\def\dropcap#1#2{%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Usage: \dropcap{1st letter}{few words}
            \settowidth{\capwd}{\dropfont#1}
            \droplinewd=\hsize%
            \advance\droplinewd by-\capwd
            %
            \advance\capwd by\dropsep
            \advance\capwd by\caphskip
            \advance\droplinewd by-\dropsep
            \advance\droplinewd by-\caphskip
            %
            \advance\droplines by1
            \vskip-\baselineskip
            \moveright \caphskip \vbox to 0pt{\vspace*{\capvskip}\hbox to \capwd{{\dropfont#1}}}\par%
            \parshape=  \droplines%
                        \ifcase\droplines%
                    \or%
                    \or%
                        \capwd \droplinewd
                        \or
                        \capwd \droplinewd \capwd \droplinewd
                        \or
                        \capwd \droplinewd \capwd \droplinewd \capwd \droplinewd
                        \or
                        \capwd \droplinewd \capwd \droplinewd \capwd \droplinewd \capwd \droplinewd
                        \fi
                        0pt \hsize
            \noindent\ignorespaces\textbf{#2}}
%
%%%%%%%%%%%%%%%%%%%% ** ARTICLE SECTIONS ** %%%%%%%%%%%%%%%%%%%
%
\def\cleardoublepage{\clearpage\if@twoside \ifodd\c@page\else
    \hbox{}\newpage\if@twocolumn\hbox{}\newpage\fi\fi\fi}
\newcommand*\chaptermark[1]{}
\newcommand*\authormark[1]{}
\setcounter{secnumdepth}{0}
\newcounter {part}
\newcounter {chapter}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\newcounter {subsubparagraph}[subsubsection]
\renewcommand \thepart {\@Roman\c@part}
\renewcommand \thechapter {\@arabic\c@chapter}
\renewcommand \thesection {\@arabic\c@section}
\renewcommand\thesubsection   {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection .\@arabic\c@subsubsection}
\renewcommand\theparagraph    {\thesubsubsection.\@arabic\c@paragraph}
\renewcommand\thesubparagraph {\theparagraph.\@arabic\c@subparagraph}
\renewcommand\thesubsubparagraph{\thesubsubsection.\@arabic\c@subsubparagraph}
\newcommand\@chapapp{\chaptername}
\def\Spellit#1{{\ifcase#1\or One\or Two\or Three\or Four\or Five\or Six\or Seven\or Eight\or Nine\or Ten\or Eleven \else\@ctrerr\fi}}

\def\Thechapter{{\ifcase\c@chapter \or 1\or 2\or 3\or 4\or 5\or 6\or 7\or 8\or 9
\or 1\hskip2.7mm0
\or 1\hskip2.7mm1\or 1\hskip2.7mm2\or 1\hskip2.7mm3\or 1\hskip2.7mm4\or 1\hskip2.7mm5\or 1\hskip2.7mm6
\or 1\hskip2.7mm7
\or 1\hskip2.7mm8\or 1\hskip2.7mm9\or 2\hskip2.7mm0\or 2\hskip2.7mm1\or 2\hskip2.7mm2\or 2\hskip2.7mm3
\or 2\hskip2.7mm4\or 2\hskip2.7mm5\or 2\hskip2.7mm6\or 2\hskip2.7mm7\or 2\hskip2.7mm8\or 2\hskip2.7mm9
\or 3\hskip2.7mm0\or 3\hskip2.7mm1\or 3\hskip2.7mm2\or 3\hskip2.7mm3\or 3\hskip2.7mm4\or 3\hskip2.7mm5
\or 3\hskip2.7mm6\or 3\hskip2.7mm7\or 3\hskip2.7mm8\or 3\hskip2.7mm9\or 4\hskip2.7mm0\or 4\hskip2.7mm1
\or 4\hskip2.7mm2\or 4\hskip2.7mm3\or 4\hskip2.7mm4\or 4\hskip2.7mm5\else\@ctrerr\fi}}


\newcommand\part{%
   \if@noskipsec \leavevmode \fi
   \par
   \addvspace{4ex}%
   \@afterindenttrue
   \secdef\@part\@spart}

\def\@part[#1]#2{%
    \ifnum \c@secnumdepth >\m@ne
      \refstepcounter{part}%
      \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
    \else
      \addcontentsline{toc}{part}{#1}%
    \fi
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont
     \ifnum \c@secnumdepth >\m@ne
       \Large\bfseries \partname~\thepart
       \par\nobreak
     \fi
     \huge \bfseries #2%
     \markboth{}{}\par}%
    \nobreak
    \vskip 3ex
    \@afterheading}
\def\@spart#1{%
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont
     \huge \bfseries #1\par}%
     \nobreak
     \vskip 3ex
     \@afterheading}
%
\def\@seccntformat#1{\csname the#1\endcsname.\hskip0.4em\relax}

\newcommand\abstracthead{\@startsection {section}{4}{1em}%
                                   {0pt \@plus 0\p@ \@minus 0\p@}%
                                   {-4pt}%
                                   {\abstractheadfont}}
\newcommand\section{\@startsection {section}{1}{\z@}%
                                   {-8pt \@plus -2\p@ \@minus -2\p@}%
                                   {4pt}%
                                   {\sectionfont}}
\newcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                   {-10pt \@plus -2\p@ \@minus -2\p@}%
                                   {2\p@}%
                                   {\subsectionfont}}
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                    {-15pt \@plus -2\p@ \@minus -2\p@}%
                                    {0.1\p@}%
                                    {\subsubsectionfont}}
\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                    {6\p@ \@plus 2\p@ \@minus 2\p@}%
                                    {-3\p@ \@plus -1\p@ \@minus -1\p@}%
                                    {\paragraphfont}}
\newcommand\subparagraph{\@startsection{subparagraph}{5}{0pt}%
                                    {0\p@ \@plus 0\p@ \@minus 0\p@}%
                                    {-.5em \@plus -1\p@ \@minus -1\p@}%
                                    {\subparagraphfont}}
\newcommand\subsubparagraph{\@startsection{subsubparagraph}{4}{\z@}%
                                    {-11pt \@plus -2\p@ \@minus -2\p@}%
                                    {-1em}%
                                    {\subsubparagraphfont}}
\newcommand\listhead{\@startsection{subparagraph}{5}{12\p@}%
                                    {0\p@ \@plus 0\p@ \@minus 0\p@}%
                                    {0.0001pt}%
                                    {\listheadfont}}

\def\periodaftersection#1{\gdef\@periodaftersection{#1}}
\periodaftersection{.}

\def\fparagraph{%
   \@startsection{paragraph}{4}{\z@}
   {-6\p@}
   {-.4em}
   {\reset@font\itshape}}

\def\@sect#1#2#3#4#5#6[#7]#8{%
    \ifnum #2>\c@secnumdepth%
        \let\@svsec\@empty%
    \else%
        \refstepcounter{#1}%
        \ifnum #2 >4%
            \edef\@svsec{\csname the#1\endcsname\hskip .25em}%
        \else%
            \edef\@svsec{\csname the#1\endcsname\hskip1em}%
        \fi%
    \fi%
    \@tempskipa #5\relax%
    \ifdim \@tempskipa>\z@%
    \@afterindentfalse
        \ifcase#2\or%
            \begingroup #6\relax
            \parindent=0pt%
            \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M{{#8}}\par}
            \endgroup
        \or
            %\@afterindenttrue
            \begingroup #6\relax
            \parindent=0pt%
            \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M{{#8}}\par}%
            \endgroup
        \or
            %\@afterindenttrue
            \begingroup #6\relax
            \parindent=0pt%
            \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M{{#8}}\par}%
            \endgroup
        \or
            \begingroup #6\relax
            \parindent=0pt%
            \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M{{#8}}\par}%
            \endgroup
        \or
            \begingroup #6\relax
            \parindent=0pt%
            \@hangfrom{\hskip #3\relax}{\interlinepenalty \@M{{#8}}\par}%
            \endgroup
        \else
            \begingroup #6\relax
            {\noindent{\hskip #3\relax\@svsec}{\interlinepenalty \@M{#8}\par}}%
            \endgroup
        \fi%
        \csname #1mark\endcsname{#7}%
        \addcontentsline{toc}{#1}{\ifnum #2>\c@secnumdepth%
                                  \else \protect\numberline{\csname the#1\endcsname}%
                                  \fi#7}%
    \else%
    \periodaftersection{}
    \def\@svsechd{%
      #6{\hskip #3\relax
      \@svsec #8\@periodaftersection}%
      \csname #1mark\endcsname{#7}%
      \addcontentsline{toc}{#1}{%
        \ifnum #2>\c@secnumdepth \else
          \protect\numberline{\csname the#1\endcsname}%
        \fi
        #7}}%
    \fi%
\@xsect{#5}}
%
\def\@ssect#1#2#3#4#5{%
  \@tempskipa #3\relax
  \ifdim \@tempskipa>\z@
  \ifcase\level
  \or
    \begingroup
      #4{%
        \@hangfrom{\hskip #1}%
          \interlinepenalty \@M \uppercase{#5}\@@par}%
    \endgroup
  \or
    %\@afterindenttrue
    \begingroup
      #4{%
        \@hangfrom{\hskip #1}%
          \interlinepenalty \@M #5\@@par}%
    \endgroup
  \or
    %\@afterindenttrue
    \begingroup
      #4{%
        \@hangfrom{\hskip #1}%
          \interlinepenalty \@M #5\@@par}%
    \endgroup
  \else
    \begingroup
      #4{%
        \@hangfrom{\hskip #1}%
          \interlinepenalty \@M #5\@@par}%
    \endgroup
  \fi
  \else
    \def\@svsechd{#4{\hskip #1\relax #5}}%
  \fi
  \@xsect{#3}}
%
\def\@startsection#1#2#3#4#5#6{%
  \if@noskipsec \leavevmode \fi%
  \par%
  \def\level{#2}%%%%%% my Change
  \@tempskipa #4\relax%
  \@afterindenttrue%
  \ifdim \@tempskipa <\z@%
    \@tempskipa -\@tempskipa \@afterindenttrue%
  \fi%
  \if@nobreak%
    \everypar{}%
  \ifnum#2=2%
  \vspace*{.5pt}%
  \fi%
  \else%
    \addpenalty\@secpenalty\addvspace\@tempskipa%
  \fi%
  \@ifstar%
    {\@ssect{#3}{#4}{#5}{#6}}%
    {\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}}
%
%%%%%%%%%%%%%%%%%%%% ** LISTS ** %%%%%%%%%%%%%%%%%%%
%
\if@twocolumn
  \setlength\leftmargini  {2em}
\else
  \setlength\leftmargini  {2.5em}
\fi
\leftmargin  \leftmargini
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\if@twocolumn
  \setlength\leftmarginv  {.5em}
  \setlength\leftmarginvi {.5em}
\else
  \setlength\leftmarginv  {1em}
  \setlength\leftmarginvi {1em}
\fi
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Alph\c@enumiv}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{(\theenumii)}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{\textbullet}
\newcommand\labelitemii{\normalfont\bfseries \textendash}
\newcommand\labelitemiii{\textasteriskcentered}
\newcommand\labelitemiv{\textperiodcentered}

\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}

\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}
\newenvironment{verse}
               {\let\\\@centercr
                \list{}{\itemsep      \z@
                        \itemindent   -1.5em%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 1.5em}%
                \item\relax}
               {\endlist}

\renewenvironment{itemize}
{\par\list
             {\labelitemi}%
             {\setlength{\labelwidth}{.6em}%
              \setlength{\labelsep}{.6em}%
              \setlength{\itemindent}{0em}%
              \setlength{\listparindent}{0em}%
              \setlength{\leftmargin}{11pt}%
              \setlength{\rightmargin}{0pt}%
              \setlength{\parsep}{0ex}%
              \setlength{\topsep}{4.5pt}%
              \setlength{\itemsep}{0pt}%
             }%
}
{\endlist}

\newenvironment{quotation}
               {\list{}{\listparindent 1.5em%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        \z@ \@plus\p@}%
                \item\relax}
               {\endlist}

\newenvironment{quote}
               {\quotefont\list{}{\rightmargin0pt\leftmargin0pt}%
                \item\relax}
               {\endlist}

\setlength\partopsep{2\p@ \@plus 1\p@ \@minus 1\p@}
\def\@listi{\leftmargin\leftmargini
            \parsep 4\p@ \@plus2\p@ \@minus\p@
            \topsep 8\p@ \@plus2\p@ \@minus4\p@
            \itemsep4\p@ \@plus2\p@ \@minus\p@}
\let\@listI\@listi
\@listi
\def\@listii {\leftmargin\leftmarginii
              \labelwidth\leftmarginii
              \advance\labelwidth-\labelsep
              \topsep    4\p@ \@plus2\p@ \@minus\p@
              \parsep    2\p@ \@plus\p@  \@minus\p@
              \itemsep   \parsep}
\def\@listiii{\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep
              \topsep    2\p@ \@plus\p@\@minus\p@
              \parsep    \z@
              \partopsep \p@ \@plus\z@ \@minus\p@
              \itemsep   \topsep}
\def\@listiv {\leftmargin\leftmarginiv
              \labelwidth\leftmarginiv
              \advance\labelwidth-\labelsep}
\def\@listv  {\leftmargin\leftmarginv
              \labelwidth\leftmarginv
              \advance\labelwidth-\labelsep}
\def\@listvi {\leftmargin\leftmarginvi
              \labelwidth\leftmarginvi
              \advance\labelwidth-\labelsep}
%
%
\newdimen\labelwidthi
\newdimen\labelwidthii
\newdimen\labelwidthiii
\newdimen\labelwidthiv
%
\newskip\topsepi
\newskip\topsepii
\newskip\topsepiii
\newskip\topsepiv
%
\newskip\itemsepi
\newskip\itemsepii
\newskip\itemsepiii
\newskip\itemsepiv
%
\newdimen\labelsepi
\newdimen\labelsepii
\newdimen\labelsepiii
\newdimen\labelsepiv
%
\newdimen\parsepi
\newdimen\parsepii
\newdimen\parsepiii
\newdimen\parsepiv
%
\newdimen\rightmargini
\newdimen\rightmarginii
\newdimen\rightmarginiii
\newdimen\rightmarginiv
%
\newdimen\itemindenti
\newdimen\itemindentii
\newdimen\itemindentiii
\newdimen\itemindentiv
%
\newdimen\labelleftmargin
%
\newif\ifdotalign\global\dotaligntrue
%
\def\myenumerate#1{%
  \ifnum \@enumdepth >\thr@@\@toodeep\else%
    \advance\@enumdepth\@ne%
    \settowidth{\csname labelwidth\romannumeral\the\@enumdepth\endcsname}{\listlabelfont#1}%
    \labelsep   \csname labelsep\romannumeral\the\@enumdepth\endcsname%
    \csname leftmargin\romannumeral\the\@enumdepth\endcsname\csname labelwidth\romannumeral\the\@enumdepth\endcsname%
    \advance\csname leftmargin\romannumeral\the\@enumdepth\endcsname\csname labelsep\romannumeral\the\@enumdepth\endcsname%
    \ifnum \@enumdepth =\@ne\advance\csname leftmargin\romannumeral\the\@enumdepth\endcsname\labelleftmargin\fi%%%%% Label leftmargin set for first level only!
    \edef\@enumctr{enum\romannumeral\the\@enumdepth}%
    \expandafter%
      \list%
        \csname label\@enumctr\endcsname%
        {\listfont%
         \topsep     \csname topsep\romannumeral\the\@enumdepth\endcsname%
         \itemsep    \csname itemsep\romannumeral\the\@enumdepth\endcsname%
         \rightmargin\csname rightmargin\romannumeral\the\@enumdepth\endcsname%
         \parsep     \csname parsep\romannumeral\the\@enumdepth\endcsname%
         \itemindent \csname itemindent\romannumeral\the\@enumdepth\endcsname
         \usecounter\@enumctr\def\makelabel##1{%
         \ifdotalign\hss\llap{\listlabelfont##1}%
         \else
         \hfill\llap{{\hbox to\csname labelwidth\romannumeral\the\@enumdepth\endcsname{{{\listlabelfont##1}}}}}%
         \fi}}%
  \fi}
\let\endmyenumerate =\endlist
%
\newenvironment{arabiclist}[1][(1)]{%
\dotaligntrue%
\topsepi6pt plus2pt minus0pt%
\topsepii6pt plus0pt minus0pt%
\itemsepi0pt plus0pt minus0pt%
\itemsepii0pt plus0pt minus0pt%
\itemsepiii0pt plus0pt minus0pt%
\labelsepi6pt%
\labelsepii6pt%
\labelsepiii6pt%
\rightmargini0pt%
\leftmarginii0pt%
\rightmarginiii0pt%
\parindent=12pt
\labelleftmargin3pt% leftmargin for labels for first level only
%\raggedright%
\def\labelenumi{\theenumi)}    \def\theenumi{\arabic{enumi}}
\def\labelenumii{(\theenumii)}  \def\theenumii{\arabic{enumii})}
\def\listfont{\normalfont}
\def\listlabelfont{\normalfont}
\myenumerate{#1}%
}{\endlist}

\newenvironment{alphalist}[1][(b)]{%
\dotalignfalse%
\topsepi5.5pt plus2pt minus0pt%
\topsepii5.5pt plus2pt minus0pt%
\itemsepi3pt plus0pt minus0pt%
\itemsepii3pt plus0pt minus0pt%
\labelsepi4pt%
\labelsepii4pt%
\rightmargini0pt
\rightmarginii0pt
\labelleftmargin\parindent% leftmargin for labels for first level only
%\raggedright%
\def\labelenumi{(\theenumi)}    \def\theenumi{\alph{enumi}}
\def\labelenumii{(\theenumii)}  \def\theenumii{\alph{enumii}}
\def\listfont{\normalfont}
\def\listlabelfont{\normalfont}
\myenumerate{#1}%
}{\endlist}

\newenvironment{romanlist}[1][(i)]{%
%\dotalignfalse%
\topsepi6pt plus2pt minus0pt%
\topsepii0pt%
\itemsepi0pt
\itemsepii0pt
\labelsepi4pt%
\labelsepii4pt%
\rightmargini0pt
\rightmarginii0pt
\labelleftmargin1em% leftmargin for labels for first level only
%\raggedright%
\def\labelenumi{\theenumi}    \def\theenumi{(\roman{enumi})}%btr650 changed by cepha%
\def\labelenumii{\theenumii}  \def\theenumii{\roman{enumii}}
\def\listfont{\normalfont}
\def\listlabelfont{\normalfont}
\myenumerate{#1}%
}{\endlist}
%
\newenvironment{greeklist}[1][(\ensuremath{\gamma})]{%
%\dotalignfalse%
\topsepi6pt plus2pt minus0pt%
\topsepii6pt plus2pt minus0pt%
\itemsepi0pt
\itemsepii0pt
\labelsepi4pt%
\labelsepii4pt%
\rightmargini0pt
\rightmarginii0pt
\labelleftmargin1em% leftmargin for labels for first level only
%\raggedright%
\def\labelenumi{\theenumi}    \def\theenumi{(\@greek\c@enumi)}
\def\labelenumii{\theenumii}  \def\theenumii{(\@greek\c@enumii)}
\def\listfont{\normalfont}
\def\listlabelfont{\normalfont}
\myenumerate{#1}%
}{\endlist}
%
\def\@greek#1{\ensuremath{\ifcase#1\or \alpha\or \beta\or \gamma\or
   \delta\or \epsilon\or \lambda\or \mu\or \sigma
   \or \zeta \else\@ctrerr\fi}}
%
\def\bulletfont{}
%
\newenvironment{bulletlist}[1][\bulletfont\textbullet]{%
\topsepi5pt plus2pt minus0pt%
\topsepii0pt%
\itemsepi7pt%
\itemsepii7pt
\labelsepi5.5pt%
\labelsepii5.5pt%
\rightmargini0pt
\rightmarginii0pt
\labelleftmargin10pt% leftmargin for labels for first level only
%\raggedright%
\def\labelenumi{\lower0pt\hbox{\bulletfont\textbullet}}%
\def\labelenumii{\lower0pt\hbox{\bulletfont\textbullet}}%
\def\listfont{\normalfont}
\def\listlabelfont{\normalfont}
\myenumerate{#1}%
}{\endlist}
%
\newenvironment{methbulletlist}[1][\bulletfont\textbullet]{%
\topsepi5pt plus2pt minus0pt%
\topsepii0pt%
\itemsepi2pt%
\itemsepii2pt
\labelsepi5.5pt%
\labelsepii5.5pt%
\rightmargini0pt
\rightmarginii0pt
\labelleftmargin9.5pt% leftmargin for labels for first level only
%\raggedright%
\def\labelenumi{\lower0pt\hbox{\bulletfont\textbullet}}%
\def\labelenumii{\lower0pt\hbox{\bulletfont\textbullet}}%
\def\listfont{\normalfont}
\def\listlabelfont{\normalfont}
\myenumerate{#1}%
}{\endlist}
%
\newenvironment{circlelist}[1][\bulletfont$\circ$]{%
\topsepi6pt plus2pt minus0pt%
\topsepii0pt%
\itemsepi0pt%
\itemsepii0pt
\labelsepi5.5pt%
\labelsepii5.5pt%
\rightmargini0pt \rightmarginii0pt
\labelleftmargin9.5pt% leftmargin for labels for first level only
%\raggedright%
\def\labelenumi{\lower0pt\hbox{\bulletfont$\circ$}}%
\def\labelenumii{\lower0pt\hbox{\bulletfont$\circ$}}%
\def\listfont{\normalfont}
\def\listlabelfont{\normalfont}
\myenumerate{#1}%
}{\endlist}
%
\newenvironment{hyphenlist}[1][\textendash]{%
\topsepi6pt plus2pt minus0pt%
\topsepii0pt plus0pt minus0pt%
\itemsepi0pt%
\itemsepii0pt
\labelsepi5.5pt%
\labelsepii5.5pt%
\rightmargini0pt
\rightmarginii0pt
\labelleftmargin1em% leftmargin for labels for first level only
%\raggedright%
\def\labelenumi{\lower0pt\hbox{\textendash}}%
\def\labelenumii{\lower0pt\hbox{\textendash}}%
\def\listfont{\normalfont}
\def\listlabelfont{\normalfont}
\myenumerate{#1}%
}{\endlist}
%
%
\newenvironment{unnumlist}[1][]{%
\topsepi6pt plus2pt minus0pt%
\topsepii0pt%
\itemsepi0pt%
\itemsepii0pt
\labelsepi0pt%
\labelsepii0pt%
\rightmargini0pt
\rightmarginii0pt
\labelleftmargin2em% leftmargin for labels for first level only
\itemindenti-1em%%%%% hangindent
%\raggedright%
\def\labelenumi{}%
\def\labelenumii{}%
\def\listfont{\normalfont}
\def\listlabelfont{\normalfont}
\myenumerate{#1}%
}{\endlist}
%
\gdef\tab#1{\begin{tabular}{#1}}
%
\newenvironment{mclist}[1]{%\par\vskip6pt plus2pt minus0pt
\parindent0pt\mclistfont
\list{}{\topsep6pt%
        \parsep0pt%
        \itemindent0pt
        \leftmargin0pt}\item[]#1}
{\end{tabular}\endlist}%
%
%%%%%%%%%%%%%%%%%%%% ** MATHS ** %%%%%%%%%%%%%%%%%%%
%
\renewcommand \theequation {\@arabic\c@equation}
\setlength\arraycolsep{1.5\p@}
\setlength\arrayrulewidth{.4\p@}

%%% To get better spacing in the subscript
%%% \ldots and

\DeclareMathSymbol{\ldotjs}{\mathord}{letters}{"3A}
\def\mymathellipsis{{\mathchoice
  {\mathinner{\mskip0.05mu\ldotjs\mskip4mu\ldotjs\mskip4mu
    \ldotjs\mskip4mu}}%
  {\mathinner{\mskip0.05mu\ldotjs\mskip4mu\ldotjs\mskip4mu
    \ldotjs\mskip4mu}}%
  {\mathinner{\mskip2.6mu\ldotjs\mskip3.3mu\ldotjs\mskip3.3mu
    \ldotjs\mskip3.3mu}}%
  {\mathinner{\mskip2.6mu\ldotjs\mskip3.3mu\ldotjs\mskip3.3mu
    \ldotjs\mskip3.3mu}}}}%
\DeclareRobustCommand{\mydots}{%
   \ifmmode\mymathellipsis\else\textellipsis\fi}%
%\AtBeginDocument{\let\ldots=\mydots}

\def\mathtight{%
\thickmuskip2mu
\thinmuskip1mu
\medmuskip=1.5mu\relax}

%%% Standard math spacing
\def\MathSpS{%
  \thinmuskip=3mu
  \medmuskip=4mu plus 2mu minus 4mu
  \thickmuskip=5mu plus 5mu\relax}
%%% With no extra declaration it is \MathSPS in effect

%%% Standard math spacing with less flexibility
\def\MathSpSS{%
  \thinmuskip=3mu
  \medmuskip=4mu plus 1.5mu minus 2mu
  \thickmuskip=5mu plus 2mu minus 1mu\relax}

%%% Tight math spacing
\def\MathSpT{%
  \thinmuskip=2mu
  \medmuskip=2.8mu plus 1.5mu minus 1.5mu
  \thickmuskip=3.7mu plus 2mu minus 1mu\relax}

%%% Very tight math spacing
\def\MathSpTT{%
  \thinmuskip=2mu
  \medmuskip=2mu plus 1.2mu minus 0.8mu
  \thickmuskip=2.6mu plus 1.5mu minus 1mu\relax}

\jot=3pt
\def\DS{\displaystyle}
\mathchardef\@m=1500   % adapted value
\mathtight
%
%%%%%%%%%%%%%%%%%%%% ** OTHER ENVIR. ** %%%%%%%%%%%%%%%%%%%
%
\newenvironment{extract}
    {\begin{list}{}
    {%
        \extractfont
        \setlength{\topsep}{6pt}
        \setlength{\leftmargin}{24pt}
        \setlength{\rightmargin}{24pt}
        \setlength{\parsep}{0pt}
        \setlength{\itemsep}{18pt}
        \setlength{\itemindent}{0pt}
        \setlength{\listparindent}{12pt}
    }\item[]}
    {\end{list}}
\def\source#1{{\par\selectfont\leftskip0pt plus1fill #1\par}}
%
\newenvironment{ack}[1][\ackname]
               {\section*{#1}%
                \ackfont\list{}{\rightmargin0pt\leftmargin0pt}%
                \item\relax}
               {\endlist}
%
%%%%%%%%%%%%%%%%%%%% ** FIGURE & TABLE ** %%%%%%%%%%%%%%%%%%%
%
\setlength\tabcolsep{6\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
%
\newcounter{figure}
\def\thefigure{\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure.}
%
\newcounter{table}
\def\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable.}
%
\newenvironment{figure}
               {\let\@makecaption\@figurecaption\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\let\@makecaption\@figurecaption\@dblfloat{figure}}
               {\end@dblfloat}
\newenvironment{intextfigure}
               {\par\addvspace{12pt plus2pt minus2pt}\begin{minipage}{\textwidth}\centering}%
               {\end{minipage}\par\addvspace{12pt plus2pt minus2pt}\@endparenv}
%
\newenvironment{table}
               {\let\@makecaption\@tablecaption\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\let\@makecaption\@tablecaption\@dblfloat{table}}
               {\end@dblfloat}
%
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{12\p@}
\setlength\belowcaptionskip{-6\p@}
%
%\long\def\@makecaption#1#2{%
%  \vskip\abovecaptionskip
%  \sbox\@tempboxa{#1: #2}%
%  \ifdim \wd\@tempboxa >\hsize
%    #1: #2\par
%  \else
%    \global \@minipagefalse
%    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
%  \fi
%  \vskip\belowcaptionskip}
%
\long\def\@makecaption#1#2{%
  \vskip0pt%\abovecaptionskip
  \sbox\@tempboxa{\noindent\figurecaptionfont {\figurenumberfont #1}\hskip0.5em#2\vphantom{pyg}}%
  \ifdim \wd\@tempboxa >\hsize
   \noindent\figurecaptionfont {\figurenumberfont #1}\hskip0.5em#2\vphantom{pyg}\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\box\@tempboxa\hfil}%
  \fi
  \vskip2pt
  \hrule height .5pt
  \vskip0pt}
%
\def\periodaftercaption#1{\gdef\@periodaftercaption{#1}}
\periodaftercaption{.}

\long\def\@figurecaption#1#2{%
\ifx\@webpdf\undefined\else%
\def\HyperRaiseLinkDefault{\figureheight}%
\advance\HyperRaiseLinkDefault\abovecaptionskip%
\fi%
  \vskip-6pt%
  \sbox\@tempboxa{\figurecaptionfont {\figurenumberfont #1\@periodaftercaption}\hskip0.5em#2}%
  \ifdim \wd\@tempboxa > \hsize
     {\figurecaptionfont {\figurenumberfont #1\@periodaftercaption}\hskip0.5em#2\vphantom{pyg}\par}
  \else
     \global \@minipagefalse
     \hb@xt@\hsize{{\unhbox\@tempboxa\hfill}}%
  \fi}
%
\graphicspath{{\artdir}{\artdir}}
%
\newdimen\figurewidth
\newdimen\figureheight%
\def\EPSfbox#1{\centering%
\settowidth{\figurewidth}{\epsfbox{#1}}%
\settoheight{\figureheight}{\epsfbox{#1}}%
\vspace*{\@fptop}%
\epsfbox{#1}
}%
%
\newdimen\dwidth
\newdimen\dheight
\newcount\ncolumns
\newcount\nrows
%
\def\showtiff#1#2#3#4{%
\dwidth=#2\dheight=#2
\ncolumns=#3\nrows=#4
\divide\dheight by \ncolumns
\multiply\dheight by \nrows
\edef\width{\number\dwidth}\edef\height{\number\dheight}\hbox{\vbox to
\dheight{\vfill\hbox to \dwidth{\special{insertimage:#1 \width \space
\height}\hfill}}}\vspace{4pt}}%
%
% TABLE
%
\newdimen\tablewd
\newdimen\tablewidth
\newdimen\captionwd
\newbox\tempbox
\newskip\tableleftskip%
\newskip\tablerightskip%
\newskip\tablecapleftskip%
\newskip\tablecaprightskip%
\newskip\tablefootleftskip%
\newskip\tablefootrightskip%
%
\newcolumntype{d}[1]{D{.}{.}{#1}}
%
\def\TCH#1{\TCHfont #1} %Table col-heads
\def\SCH#1{\SCHfont #1} %Table col-heads
%
\def\periodaftercaption#1{\gdef\@periodaftercaption{#1}}%S%
\periodaftercaption{}%S%
\def\@tablecaption#1#2{%
  \vspace*{\@fptop}%
  \settowidth{\global\captionwd}{\hbox{\tablecaptionfont{\tablenumberfont #1.}\hskip.5em#2}}%
  {\tablecaptionfont{\tablenumberfont #1\@periodaftercaption}\hskip.5em#2\par}
  \vskip\belowcaptionskip}
%
\long\def\tableparts#1#2#3{%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #1-caption #2-body #3-footnote
    \setbox\tempbox\vbox{#1}                          %% to calculate width caption
    \addtocounter{table}{-1}                          %%
    \settowidth{\global\tablewd}{\tablebodyfont#2}    %% \captionwd - caption of width
    \tablewidth\hsize                                 %% \tablewd   - table width
    \advance\tablewidth-\tablewd
    \divide\tablewidth 2
    \global\tableleftskip\tablewidth
    \global\tablerightskip\tablewidth
    \ifdim\captionwd>\hsize
        \global\tablecapleftskip\tablewidth
        \global\tablecaprightskip\tablewidth %plus2fill
    \else
        \global\tablecapleftskip\tablewidth  %plus2fill
        \global\tablecaprightskip\tablewidth plus2fill
    \fi
    \global\tablefootleftskip\tablewidth
    \global\tablefootrightskip\tablewidth %plus2fill
    {#1 \par}%
    {\tablebodyfont #2\par}%
    {\tablefootfont #3\par}%
\vskip6pt}
%B%\long\def\tableparts#1#2#3{%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #1-caption #2-body #3-footnote
%B%    \setbox\tempbox\vbox{#1}                          %% to calculate width caption
%B%    \addtocounter{table}{-1}                          %%
%B%    \settowidth{\global\tablewd}{\tablebodyfont#2}    %% \captionwd - caption of width
%B%    \tablewidth\hsize                                 %% \tablewd   - table width
%B%    \advance\tablewidth-\tablewd
%B%    \divide\tablewidth 2
%B%    \global\tableleftskip\tablewidth
%B%    \global\tablerightskip\tablewidth
%B%    \ifdim\captionwd>\hsize
%B%        \global\tablecapleftskip\tablewidth
%B%        \global\tablecaprightskip\tablewidth plus2fill
%B%    \else
%B%        \global\tablecapleftskip\tablewidth  %plus2fill
%B%        \global\tablecaprightskip\tablewidth plus2fill
%B%    \fi
%B%    \global\tablefootleftskip\tablewidth
%B%    \global\tablefootrightskip\tablewidth %plus2fill
%B%    {#1 \par}%
%B%    {\tablebodyfont #2\par}%
%B%    \vspace{3pt}
%B%    {\tablefootfont #3\par}%
%B%}
%
%
\long\def\mytableparts#1#2#3{%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #1-caption #2-body #3-footnote
    \setbox\tempbox\vbox{#1}                          %% to calculate width caption
    \addtocounter{table}{-1}                          %%
    \settowidth{\global\tablewd}{\tablebodyfont#2}    %% \captionwd - caption of width
    \tablewidth\hsize                                 %% \tablewd   - table width
    \advance\tablewidth-\tablewd
    \divide\tablewidth 2
    \global\tableleftskip\tablewidth
    \global\tablerightskip\tablewidth
    \ifodd\thepage\else\advance\tableleftskip4pt\fi
    \ifdim\captionwd>\hsize
        \global\tablecapleftskip\tablewidth
        \global\tablecaprightskip\tablewidth %plus2fill
    \else
        \global\tablecapleftskip\tablewidth  %plus2fill
        \global\tablecaprightskip\tablewidth plus2fill
    \fi
    \ifodd\thepage\else\advance\tablecapleftskip4pt\fi
    \global\tablefootleftskip\tablewidth
    \global\tablefootrightskip\tablewidth %plus2fill
    \ifodd\thepage\else\advance\tablefootleftskip4pt\advance\tablefootrightskip4pt\fi
    {#1 \par}%
    {\tablebodyfont #2\par}%
    {\tablefootfont #3\par}%
}
%
% Tablenotes
\long\def\tablenote#1{\par #1\par}
\long\def\tablesource#1{\par\textit{Source}:\ #1\par}
\long\def\tablefoot#1{\par#1\par}
%
%
\newif\if@AProtate\@AProtatefalse
\def\AProtate{\global\@AProtatetrue}
\def\endAProtate{\global\@AProtatefalse}
\newdimen\rotdimen
\def\rotstart#1{\special{ps: gsave currentpoint currentpoint translate
        #1 neg exch neg exch translate}}
\def\rotfinish{\special{ps: currentpoint grestore moveto}}
\def\rotl#1{\rotdimen=\ht#1\advance\rotdimen by \dp#1
        \hbox to \rotdimen{\vbox to \wd#1{\vskip \wd#1
        \rotstart{-270 rotate}\box #1\vss}\hss}\rotfinish}
\def\rotr#1{\rotdimen=\ht#1\advance\rotdimen by \dp#1
        \hbox to \rotdimen{\vbox to \wd#1{\vskip \wd#1
        \rotstart{90 rotate}\box #1\vss}\hss}\rotfinish}
\def\Rotate#1{\setbox4=\hbox to\textheight{\vbox{\hsize\textheight
\begin{minipage}[t]{\textheight}
\begin{center}#1\end{center}
\end{minipage}}}\centerline{\rotl{4}}}
%
%
\def\tfmark#1{\if#1**\else\textsuperscript{\smash{#1}}\fi}
%
\def\myhline#1{%
  \noalign{\ifnum0=`}\fi\hrule \@height #1 \futurelet
   \reserved@a\@xhline}%
%
%\def\toprule{\myhline{1pt}\vrule height12pt depth0pt width0pt}%
%\def\colrule{\vrule height0pt depth5pt width0pt\\\myhline{1pt}\vrule height12pt depth0pt width0pt}%
%\def\botrule{\vrule height0pt depth5pt width0pt\\\myhline{1pt}\omit {\vspace*{4pt}}}%
%

\newcommand{\up}{\HD{10}{0}}
\newcommand{\down}{\HD{0}{5}}
\newcommand{\Ndown}{\HD{0}{5}}

\def\toprule{\noalign{\ifnum0=`}\fi
  \hrule \@height 1\p@
  \hrule \@height 2.5\p@ \@width \z@
  \futurelet \@tempa\@xhline\omit}

\def\colrule{\noalign{\ifnum0=`}\fi
  \hrule \@height .5\p@ \@width \z@
  \hrule \@height .5\p@
  \hrule \@height 2.5\p@ \@width\z@
  \futurelet \@tempa\@xhline\omit}

\def\botrule{\noalign{\ifnum0=`}\fi
  \hrule \@height -.5\p@ \@width \z@
  \hrule \@height 1\p@
  \hrule \@height -4\p@ \@width\z@
  \futurelet \@tempa\@xhline\omit }

%\def\crules#1{\vrule height0pt depth5pt width0pt\\#1\vrule height11pt depth0pt width0pt}

\def\crules#1{\noalign{\hrule \@height 5.5\p@ \@width \z@}#1\noalign{\hrule \@height 5.5\p@ \@width \z@}}
%
%\def\@cline#1-#2\@nil{%
%  \omit
%  \@multicnt#1%
%  \advance\@multispan\m@ne
%  \ifnum\@multicnt=\@ne\@firstofone{&\omit}\fi
%  \@multicnt#2%
%  \advance\@multicnt-#1%
%  \advance\@multispan\@ne
%  \leaders\hrule\@height\arrayrulewidth\hfill
%  \cr
%  \noalign{\vskip-\arrayrulewidth}}
%
%%%%%%%%%%%%%%%%%%%% ** FLOAT PLACEMENT ** %%%%%%%%%%%%%%%%%%%
%
% Float placement parameters.
%
% The total number of floats that can be allowed on a page.
\setcounter{totalnumber}{9}
% The maximum number of floats at the top and bottom of a page.
\setcounter{topnumber}{9}
\setcounter{bottomnumber}{9}
% The maximum part of the top or bottom of a text page that can be
% occupied by floats. This is set so that at least four lines of text
% fit on the page.
\renewcommand\topfraction{.9}
\renewcommand\bottomfraction{.9}
% The minimum amount of a text page that must be occupied by text.
% This should accomodate four lines of text.
\renewcommand\textfraction{.1}
% The minimum amount of a float page that must be occupied by floats.
\renewcommand\floatpagefraction{.85}%{.91}
% The same parameters repeated for double column output
\setcounter{dbltopnumber}{2}
\renewcommand\dbltopfraction{.9}
\renewcommand\dblfloatpagefraction{.9}
% Space between floats
\setlength\floatsep{18\p@ \@plus 3\p@ \@minus \z@}
% Space between floats and text
\setlength\textfloatsep{11\p@ \@plus 3\p@ \@minus 3\p@}
% Space above and below an inline figure
\setlength\intextsep   {6\p@ \@plus \z@ \@minus 2\p@}
% For double column floats
\setlength\dblfloatsep    {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dbltextfloatsep{17\p@ \@plus 6\p@ \@minus 6\p@}

% Space left at top, bottom and inbetween floats on a float page.
\setlength\@fptop{3.5\p@}% \@plus 1fil         % no space above float page figures
\setlength\@fpsep{12\p@ \@plus 2fil}
\setlength\@fpbot{0\p@ \@plus 1fil}

% The same for double column
\setlength\@dblfptop{0\p@}
\setlength\@dblfpsep{12\p@ \@plus 2fil}
\setlength\@dblfpbot{0\p@ \@plus 1fil}
%
%%%%%%%%%%%%%%%%%%%% ** TOC ** %%%%%%%%%%%%%%%%%%%
%
\newcommand\@pnumwidth{1.5em}
\newcommand\@tocrmarg{0em}
\newcommand\@dotsep{2}
\setcounter{tocdepth}{1}
\def\dotfill{\leaders\hbox to.3em{\hss .\hss}\hskip\z@ plus  1fill}%
%
\newcommand\tableofcontents{%
    \section*{\contentsname
        \@mkboth{%
           \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
    \@starttoc{toc}%
    }
\newcommand*\l@part[3]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty\@secpenalty
    \addvspace{.5em \@plus\p@}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \@hangfrom{\hbox to 10.5pc{\bfseries #1\hfill}}\textbf{#2} \hfil \hb@xt@\@pnumwidth{\hss #3}}\par
       \nobreak
       \if@compatibility
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}%
      \fi
    \endgroup
  \fi}
\newcommand*\l@subpart[3]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty\@secpenalty
    \addvspace{.5em \@plus\p@}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \@hangfrom{\hbox to 10.5pc{\bfseries #1\hfill}}\textbf{#2} \hfil \hb@xt@\@pnumwidth{\hss #3}}\par
       \nobreak
       \if@compatibility
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}%
      \fi
    \endgroup
  \fi}
\newcommand*\l@section[2]{%
  \ifnum \c@tocdepth >\z@
    \addpenalty\@secpenalty
    \addvspace{1.0em \@plus\p@}%
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
    \endgroup
  \fi}
%
\newcommand*\l@article[2]{%
\renewcommand\@pnumwidth{5pc}
  \ifnum \c@tocdepth >\z@
    \addpenalty\@secpenalty
    \addvspace{.6em \@plus\p@}%
    \setlength\@tempdima{10.5pc}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \fontsize{9}{11}\selectfont
      \advance\leftskip\@tempdima
      %\hskip -\leftskip
      {\bfseries #1}\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
    \endgroup
  \fi}
%
\newcommand*\l@author[2]{%
  \ifnum \c@tocdepth >\z@
    \addpenalty\@secpenalty
    \addvspace{0em \@plus\p@}%
    \setlength\@tempdima{10.5pc}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \fontsize{8}{9}\selectfont
      \advance\leftskip\@tempdima
      %\hskip -\leftskip
      {#1}\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
    \endgroup
  \fi}
%

\def\@dottedtocline#1#2#3#4#5{%
  \ifnum #1>\c@tocdepth \else
    \vskip \z@ \@plus.2\p@
    {\leftskip #2\relax \rightskip \@tocrmarg \parfillskip -\rightskip
     \parindent #2\relax\@afterindenttrue
     \interlinepenalty\@M
     \leavevmode
     \@tempdima #3\relax
     \advance\leftskip \@tempdima \null\nobreak\hskip -\leftskip
     {#4}\nobreak
     \leaders\hbox{$\m@th
        \mkern \@dotsep mu\hbox{.}\mkern \@dotsep
        mu$}\hfill
     \nobreak
     \hb@xt@\@pnumwidth{\hfil\normalfont \normalcolor #5}%
     \par}%
  \fi}
%
\newcommand*\l@subsection{\@dottedtocline{2}{1.5em}{2.3em}}
\newcommand*\l@subsubsection{\@dottedtocline{3}{3.8em}{3.2em}}
\newcommand*\l@paragraph{\@dottedtocline{4}{7.0em}{4.1em}}
\newcommand*\l@subparagraph{\@dottedtocline{5}{10em}{5em}}
\newcommand\listoffigures{%
    \section*{\listfigurename
      \@mkboth{\MakeUppercase\listfigurename}%
              {\MakeUppercase\listfigurename}}%
    \@starttoc{lof}%
    }
\newcommand*\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\listoftables{%
    \section*{\listtablename
      \@mkboth{%
          \MakeUppercase\listtablename}{\MakeUppercase\listtablename}}%
    \@starttoc{lot}%
    }
\let\l@table\l@figure
%
%%%%%%%%%%%%%%%%%%%% ** FOOTNOTES ** %%%%%%%%%%%%%%%%%%%
%
\setlength\footnotesep{6.65\p@}
\setlength{\skip\footins}{11\p@ \@plus 4\p@ \@minus 0\p@}
\skip\@mpfootins = \skip\footins
%
\renewcommand\footnoterule{%
  \kern-2.3\p@
  \hrule\@width\columnwidth \@height .5\p@
  \kern2.3\p@}
%
\def\footnum#1{\gdef\@footnum{#1}}
\footnum{\mbox{\fontsize{5.31}{5}\selectfont\arabic{footnote}}}
%
\newdimen\footnumskip
\def\@makefnmark{$^{\mbox{\fontsize{5.98}{6}\selectfont\arabic{footnote}}}$}
\long\def\@makefntext#1{\hspace*{\footnumskip}{\textsuperscript{\@footnum}}#1\footnumskip=0pt}
%
%%%%%%%%%%%%%%%%%%%% ** APPENDIX ** %%%%%%%%%%%%%%%%%%%
%
\newcommand\appendix{\par
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \setcounter{table}{0}%
  \setcounter{figure}{0}%
  \setcounter{equation}{0}%
  \gdef\thesection{\Alph{section}.1}%
  \def\thefigure{\Alph{section}\@arabic\c@figure}%
  \def\thetable{\Alph{section}\@arabic\c@table}%
  \def\theequation {\Alph{section}\@arabic\c@equation}}
%
%%%%%%%%%%%%%%%%%%%% ** BIBLIOGRAPHY ** %%%%%%%%%%%%%%%%%%%
%
\newdimen\bibindent
\setlength\bibindent{11pt}
\def\@mybiblabel#1{\llap{\hbox to \leftmargin{\hfill#1\hskip\labelsep}\hskip-\labelsep}}
%
\newdimen\maxbibwd
\gdef\maxbib#1{%
\global\settowidth\maxbibwd{#1}
\leftskip-\maxbibwd}
\def\resetmaxbib{\leftskip0pt}
%
\newenvironment{thebibliography}[1]
     {\bibsection%
     \gdef\@bibmax{#1}%
     \addcontentsline{toc}{section}{\protect\numberline{\relax}\refname}
        \list{\ifx\@bibmax\empty\@biblabel{\@arabic\c@enumiv}\else\@mybiblabel{\@biblabel{\@arabic\c@enumiv}}\fi}%
           {\bibfont%
            \ifx\@bibmax\empty
                \settowidth\labelwidth{\@biblabel{#1}}%
                \leftmargin\labelwidth
                \advance\leftmargin\labelsep
            \else
                \itemindent=0pt
                \labelsep=6pt
                \settowidth\bibindent{[\@bibmax]}
                \advance\bibindent\labelsep
                \leftmargin\bibindent
            \fi
            \itemsep=0pt \@plus 0\p@
            \parsep=0pt
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
    \endlist
    }
%
\newcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
\newcommand{\bysame}{\leavevmode\hbox to5em{\hrulefill}\thinspace}
%
%%%%%%%%%%%%%%%% ** CITATION STYLE ** %%%%%%%%%%%%%%%%
\if@openbib
\DeclareRobustCommand\cite{%
  \@ifnextchar [{\@tempswatrue\@citex}{\@tempswafalse\@citex[]}}
\def\@citex[#1]#2{%
  \let\@citea\@empty
  \@cite{\@for\@citeb:=#2\do
    {\@citea\def\@citea{;\penalty\@m\ }%
     \edef\@citeb{\expandafter\@firstofone\@citeb\@empty}%
     \if@filesw\immediate\write\@auxout{\string\citation{\@citeb}}\fi
     \@ifundefined{b@\@citeb}{\mbox{\reset@font\bfseries ?}%
       \G@refundefinedtrue
       \@latex@warning
         {Citation `\@citeb' on page \thepage \space undefined}}%
       {{\csname b@\@citeb\endcsname}}}}{#1}}
\def\@cite#1#2{{#1\if@tempswa ,\ #2\fi}}
\def\@biblabel#1{}
\else
\usepackage[authoryear,sort]{cpl-natbib}
\let\mycite\citet
\let\cite\citep
\setlength{\bibhang}{13pt}
\bibsep=0pt
\fi
\def\bibfont{\fontsize{9}{18}\selectfont}
%
%B%\if@openbib
%B%\DeclareRobustCommand\cite{%
%B%  \@ifnextchar [{\@tempswatrue\@citex}{\@tempswafalse\@citex[]}}
%B%\def\@citex[#1]#2{%
%B%  \let\@citea\@empty
%B%  \@cite{\@for\@citeb:=#2\do
%B%    {\@citea\def\@citea{;\penalty\@m\ }%
%B%     \edef\@citeb{\expandafter\@firstofone\@citeb\@empty}%
%B%     \if@filesw\immediate\write\@auxout{\string\citation{\@citeb}}\fi
%B%     \@ifundefined{b@\@citeb}{\mbox{\reset@font\bfseries ?}%
%B%       \G@refundefinedtrue
%B%       \@latex@warning
%B%         {Citation `\@citeb' on page \thepage \space undefined}}%
%B%       {{\csname b@\@citeb\endcsname}}}}{#1}}
%B%\def\@cite#1#2{{#1\if@tempswa ,\ #2\fi}}
%B%
%B%\def\@biblabel#1{}
%B%\else \if@mynatbib
%B%\usepackage{natbib}
%B%\let\mycite\citet
%B%\let\cite\citep
%B%\setlength{\bibhang}{10pt}
%B%\bibsep=0pt
%B%\def\bibfont{\fontsize{8}{9}\selectfont}
%B%\else
%B%\usepackage[authoryear,sort]{cpl-natbib}
%B%\let\mycite\citet
%B%\let\cite\citep
%B%\setlength{\bibhang}{10pt}
%B%\bibsep=0pt
%B%\fi\fi
%
\def\bibfont{\fontsize{9}{18}\selectfont}
\def\bibsection{\section*{\refname}%\markboth{\refname}{\refname}
\vspace*{1pt}}
\ifx\@webpdf\undefined\else%
\usepackage{hypernat}%% used with natbib pdfstartview={FitH -32768},
\fi
%
\ifx\@webpdf\undefined
%\def\url#1{#1}%
\def\href#1#2{#1}%
\else%raiselinks,
\usepackage[hyperfootnotes, dvipsone, breaklinks=true, colorlinks=true, urlcolor=black, linkcolor=blue, citecolor=blue,linkbordercolor=0 0 0,filebordercolor=0 0 0,citebordercolor=0 0 0,
urlbordercolor=0 0 0,bookmarks,bookmarksopen=true,bookmarksnumbered=true,pdfpagelayout={OneColumn},pdfstartview={FitH -32768}]{cpl-hyperref}%
\fi%
%
%%%%%%%%%%%%%%%%%%%% ** INDEX ** %%%%%%%%%%%%%%%%%%%
%
\newenvironment{theindex}
               {%\if@twocolumn
                %  \@restonecolfalse
                %\else
                %  \@restonecoltrue
                %\fi
                \columnseprule \z@
                \columnsep 12\p@
                %\twocolumn[\section*{\indexname}]%
                \@mkboth{\MakeUppercase\indexname}%
                        {\MakeUppercase\indexname}%
                \thispagestyle{plain}\parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \let\item\@idxitem
                \linewidth42pc\hsize42pc%
                \multicolsep = 0pt plus 0pt minus 0pt%
                \begin{multicols}{3}[\vspace*{2pc}]%
                }
               {\end{multicols}%\if@restonecol\onecolumn\else\clearpage\fi
               }
\newcommand\@idxitem{\par\hangindent 18\p@}
%\newcommand\subitem{\par\hangindent 25\p@ \hspace*{8\p@}}
\newcommand\subitem{\@idxitem \hspace*{9\p@}}
\newcommand\subsubitem{\@idxitem \hspace*{18\p@}}
\newcommand\indexspace{\par \vskip 16\p@ \@plus5\p@ \@minus3\p@\relax}
%
%%%%%%%%%%%%%%%%%%%% ** OTHER MACAROS ** %%%%%%%%%%%%%%%%%%%
%
\newif\iffirstpagegiven \firstpagegivenfalse
\newif\iflastpagegiven   \lastpagegivenfalse

% use \firstpage and \lastpage to set the page range explicitly
\newcommand\firstpage[1]{%
  \gdef\@firstpage{#1}
  \global\firstpagegiventrue%
}
\def\@firstpage{1}

\newcommand\lastpage[1]{%
  \xdef\@lastpage{#1}%
  \global\lastpagegiventrue}
\def\@lastpage{0}

\def\setlastpage{\iflastpagegiven\else
  \edef\@tempa{@lastpage}%
  \expandafter
  \ifx \csname \@tempa \endcsname \relax
    \gdef\@lastpage{0}%
  \else
    \xdef\@lastpage{\@nameuse{@lastpage}}%
  \fi
  \fi}

\def\writelastpage{%
  \iflastpagegiven
  \else
    \immediate\write\@auxout%
    {\string\global\string\@namedef{@lastpage}{\thepage}}
  \fi}

\def\setfirstpage{%
  \iffirstpagegiven
  \else
    \xdef\@firstpage{\thepage}
  \fi}
%
\AtBeginDocument{\setfirstpage\setlastpage}
\AtEndDocument{\writelastpage%
%     \ifx\@history\undefined\else%
%     \vskip.5\baselineskip%
%     \vskip-.3pt%
%       {\noindent\@history\par}%
%     \fi
}
%
\newcommand\ourcurrentdirectory{}
\def\filecontents{\chaptername\space\thechapter}
\AtBeginDocument{\filename@parse{\jobname}
  \gdef\filename{\MakeLowercase{\filename@base}%
    \if\filename@ext\relax.tex%
    \else\MakeLowercase{\filename@ext}\fi}}
% calculate the time of day
\newcount\hours
\newcount\minutes
\hours=\time
\divide\hours by 60
\minutes=-\hours
\multiply\minutes by 60
\advance\minutes by \time
\def\timenow{\two@digits{\the\hours}:\two@digits{\the\minutes}}
\font\slugfont=cmr10 at 6pt
%
\IfFileExists{crop.cfg}{}{%
  \renewcommand*\CROP@@info{\slugfont%
    \lower 19\p@
    \hbox to \paperwidth{%
      \hskip 3mm
      %\scriptsize       % <--- NFSS font commands don't work here -- why?
      [\timenow\enskip\the\day/\the\month/\the\year\enskip MOLBEV-\ourcurrentdirectory\filename]
      %\hfill
      %Job No:\space\Jobid
      \hfill%\@shortjname:
      \space%\@jname
      \hfill Page:\space\thepage\quad\@firstpage--\@lastpage%
      \hskip 3mm
      }}
  \renewcommand*\CROP@@ulc{%
    \begin{picture}(0,0)\unitlength 1mm\thinlines
      \put(-13,0){\line(1,0){9}}
      \put(0,13){\line(0,-1){9}}
      \put(10,6){\copyeditor}
      \put(100,6){\msc}
    \end{picture}}
  \renewcommand*\CROP@@urc{%
    \begin{picture}(0,0)\unitlength 1mm\thinlines
      \put(13,0){\line(-1,0){9}}
      \put(0,13){\line(0,-1){9}}
    \end{picture}}
  \renewcommand*\CROP@@llc{%
    \begin{picture}(0,0)\unitlength 1mm\thinlines
      \put(-13,0){\line(1,0){9}}
      \put(0,-13){\line(0,1){9}}
    \end{picture}}
  \renewcommand*\CROP@@lrc{%
    \begin{picture}(0,0)\unitlength 1mm\thinlines
      \put(13,0){\line(-1,0){9}}
      \put(0,-13){\line(0,1){9}}
    \end{picture}}
  \cropdef\CROP@@ulc\CROP@@urc\CROP@@llc\CROP@@lrc{cam}
  \renewcommand*\CROP@kernel{\color@setgroup
    \vbox to\z@{\vskip-1in
      \hb@xt@\z@{\hskip-1in
        \CROP@every
        \vbox to\paperheight{%
          \hb@xt@\paperwidth{%
            \CROP@ulc\CROP@uedge\CROP@urc}%
          \CROP@ledge
          \hb@xt@\paperwidth{%
            \setbox\z@\hbox{\CROP@info}%
            \ht\z@\z@ \dp\z@\z@ \wd\z@\z@ \box\z@
            \CROP@llc\hfil\CROP@lrc}}%
        \hss}%
      \vss}\color@endgroup}
  }

\if@cropmarkson \crop[cam] \fi
%

%%%%%%%%%b-KOO
\long\def\koo#1#2#3{\vbox to0pt{\hsize0pt\kern #2\hbox to0pt{\kern#1{#3}\hss}\vss}}
\long\def\zeroCC#1{\vbox to0pt{\vss\hbox to0pt{\hss #1\hss}\vss}}
\long\def\zeroLC#1{\vbox to0pt{\vss\hbox to0pt{\hss #1\hss}}}
\long\def\zeroUC#1{\vbox to0pt{\hbox to0pt{\hss #1\hss}\vss}}
\long\def\zeroLR#1{\vbox to0pt{\vss\hbox to0pt{\hss #1}}}
\long\def\zeroCR#1{\vbox to0pt{\vss\hbox to0pt{\hss #1}\vss}}
\long\def\zeroUR#1{\vbox to0pt{\hbox to0pt{\hss #1}\vss}}
\long\def\zeroLL#1{\vbox to0pt{\vss\hbox to0pt{#1\hss}}}
\long\def\zeroUL#1{\vbox to0pt{\hbox to0pt{#1\hss}\vss}}
\long\def\zeroCL#1{\vbox to0pt{\vss\hbox to0pt{#1\hss}\vss}}
\def\VRHDW#1#2#3{\vrule height #1 depth #2 width #3}
\def\HRHDW#1#2#3{\hrule height #1 depth #2 width #3}
\def\HD#1#2{\vrule height #1pt depth #2pt width 0pt\relax}
%
\newdimen\grid@sidemargin
\newdimen\grid@wd
\newdimen\grid@ht
\newdimen\grid@linewd\grid@linewd=0.1pt
\newdimen\grid@topofline
\grid@topofline\headheight%
\advance\grid@topofline\headsep%
\advance\grid@topofline.35\topskip%
%
\def\oddevenswitch{\ifodd\c@page \grid@sidemargin=\oddsidemargin\else\grid@sidemargin=\evensidemargin\fi}%
%
\def\DRAWGRID{%
\def\GRIDRULES{\oddevenswitch
\grid@wd=\textwidth \advance\grid@wd by 2\grid@linewd
\grid@ht=\textheight
\advance\grid@ht\headheight
\advance\grid@ht\headsep
\advance\grid@ht by 2\grid@linewd
\advance\grid@ht by -2.1pt
\hbox to0pt{%
   \vbox to0bp{\kern-\headheight%
   \hbox to\textwidth{%
     \koo{-\grid@linewd}{0pt}{\zeroLL{\VRHDW{\grid@linewd}{0pt}{\grid@wd}}}%
     \koo{-\grid@linewd}{\headheight}{\zeroLL{\VRHDW{\grid@linewd}{0pt}{\grid@wd}}}%
     \koo{-\grid@linewd}{-\grid@linewd}{\zeroUL{\VRHDW{\grid@ht}{0pt}{\grid@linewd}}}%
     \koo{\textwidth}{-\grid@linewd}{\zeroUL{\VRHDW{\grid@ht}{0pt}{\grid@linewd}}}%
     \koo{-\grid@linewd}{\grid@ht}{\zeroUL{\VRHDW{\grid@linewd}{0pt}{\grid@wd}}}%
     \koo{-\grid@linewd}{\grid@topofline}{\zeroUL{\VRHDW{\grid@linewd}{0pt}{\grid@wd}}}%
     \hss}\vss}\hss}}%

\def\FIRSTLINEGRID{%
\ifdim\firstlinefixedpos>0pt%
\oddevenswitch%
\hbox to0pt{%
   \vbox to0bp{\kern240pt%
\hbox to\textwidth{\koo{-\grid@linewd}{10\baselineskip}{\zeroUL{\VRHDW{\grid@linewd}{0pt}{\grid@wd}}}}}}%
\fi}%
}%
%
\def\GRIDRULES{}
\def\FIRSTLINEGRID{}
%
%%%%%%%%%e-KOO
%
%%%%%%%%%%%%%%%%%%% ** LETTER SPACING & TRACKING %%%%%%%%%%%%%%%%%%
%
%
\sodef\textso{}{-.015em}{.3em\@plus.08em\@minus.06em}%
        {.55em\@plus.275em\@minus.183em}%
%
\def\dolist{\afterassignment\dodolist\let\next= }
\def\dodolist{\ifx\next\endlist \let\next\relax
   \else \\\let\next\dolist \fi
   \next}
%
\def\track#1#2{%
\def\\{\expandafter\if\space\next\setbox0=\hbox{\ }%
  \else \setbox0=\hbox{\next}\fi%
     \box0\kern#1}%
\leavevmode\hbox{\dolist#2\endlist\unkern}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SMALL CAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newdimen\currfontsize
\def\makesc#1#2{%
\currfontsize=1em
\def\scsize{\fontsize{#1\currfontsize}{1em}\selectfont}%
\def\\{\expandafter\if\space\next\gdef\xxx{\ }%
  \else\gdef\eachchar{%
  \ifx\next a{\scsize A}\else%
  \ifx\next b{\scsize B}\else%
  \ifx\next c{\scsize C}\else%
  \ifx\next d{\scsize D}\else%
  \ifx\next e{\scsize E}\else%
  \ifx\next f{\scsize F}\else%
  \ifx\next g{\scsize G}\else%
  \ifx\next h{\scsize H}\else%
  \ifx\next i{\scsize I}\else%
  \ifx\next j{\scsize J}\else%
  \ifx\next k{\scsize K}\else%
  \ifx\next l{\scsize L}\else%
  \ifx\next m{\scsize M}\else%
  \ifx\next n{\scsize N}\else%
  \ifx\next o{\scsize O}\else%
  \ifx\next p{\scsize P}\else%
  \ifx\next q{\scsize Q}\else%
  \ifx\next r{\scsize R}\else%
  \ifx\next s{\scsize S}\else%
  \ifx\next t{\scsize T}\else%
  \ifx\next u{\scsize U}\else%
  \ifx\next v{\scsize V}\else%
  \ifx\next w{\scsize W}\else%
  \ifx\next x{\scsize X}\else%
  \ifx\next y{\scsize Y}\else%
  \ifx\next z{\scsize Z}\else%
  \next\kern-.25pt\unskip\fi\fi\fi\fi\fi\fi\fi\fi\fi%
  \fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi\fi%
  }\fi\eachchar}%
\leavevmode\hbox{\dolist#2\endlist\unkern}}
%
%%%%%%%%%%%%%%%%%%%% ** Author's COMMANDS ** %%%%%%%%%%%%%%%%%%%
%
%%%%for smaller size text
\newenvironment{methods}{%
  \begingroup
  \let\bulletlist\methbulletlist%
  \let\endbulletlist\endmethbulletlist%
%\def\section{%
%  \@startsection{section}{1}{\z@}
%  {-11\p@ plus -3\p@}{4\p@}
%  {\reset@font\raggedright\helveticabold\fontsize{10}{12}\selectfont\MakeUppercase}}
% \def\subsection{%
%  \@startsection{subsection}{2}{\z@}
%  {-5\p@ plus -2\p@}{4\p@}
%  {\reset@font\raggedright\mathversion{bold}\fontseries{b}\fontsize{10}{12}\selectfont}}
% \def\subsubsection{%
%  \@startsection{subsubsection}{3}{\z@}
%  {-6\p@ plus -1\p@}{-1em}
%  {\reset@font\normalfont\normalsize\itshape}}
\footnotesize
  \par}
{\par\endgroup}%\bigskip\@afterheading\@afterindenttrue
%
\setlength{\jot}{5\p@}
\mathchardef\@m=1500   % adapted value
%
\language=2
\hyphenation{Figure Table Figures Tables}

\def\sbond{\ensuremath{\raise.25ex\hbox{${-}\!\!\!\!{-}$}}\kern -.9pt}
\def\dbond{\ensuremath{\raise.25ex\hbox{=$\!$=}}}
\def\tbond{\ensuremath{\raise.20ex\hbox{${\equiv}\!\!\!{\equiv}$}}}

%
%%%%%%%%%%%%%%%%%%%% ** USER COMMANDS ** %%%%%%%%%%%%%%%%%%%
%
\def\x{\end{document}}
\def\n{\nobreak}
\def\h{\hbox}
\def\vs#1{\vadjust{\vspace{#1pt}}}
\def\tra#1{\spaceskip 0.#1em plus 0.15em minus 0.09em\relax}%
\def\rull#1{\vspace*{#1pt}\noindent$\hrulefill$}
\def\nx{\hspace*{8pt}}
\def\SPT{\spaceskip 0.3em plus 0.16em minus 0.09em\relax}
\def\vvp{\vadjust{\vfill\pagebreak}}
\def\vvsp#1{\vadjust{\vspace{#1pt}\pagebreak}}
\def\vvSP#1{\HD{0}{#1}\pagebreak}
%
\twocolumn
\pagestyle{myheadings}
\pagenumbering{arabic}
\sloppy
\flushbottom
%
\usepackage{amsthm}
% Theorems

\def\th@plain{%
%%  \let\thm@indent\noindent % no indent
\thm@headfont{\scshape}% heading font is bold%\quad
\thm@notefont{\upshape\mdseries}% same as heading font
\thm@headpunct{.}% no period after heading
\thm@headsep 5\p@ plus\p@ minus\p@\relax%
%%  \let\thm@swap\@gobble
%%  \thm@preskip\topsep
%%  \thm@postskip\theorempreskipamount
\@bodyfont % body font
}

\def\bodyfont#1{\gdef\@bodyfont{#1}}
\bodyfont{\itshape}

\def\thmhead@plain#1#2#3{%
  \thmname{#1}\thmnumber{\@ifnotempty{#1}{ }\@upn{#2}}%
  \thmnote{ {\the\thm@notefont(#3)}}}
\def\swappedhead#1#2#3{%
  \thmnumber{#2}%
  \thmname{\@ifnotempty{#2}{~}#1}%
  \thmnote{ {\the\thm@notefont(#3)}}}
\newcounter{defn}
\def\thedefn{\arabic{defn}}
\def\periodafterdefn#1{\gdef\@periodafterdefn{#1}}
\periodafterdefn{}
\def\enddefn#1{\gdef\@enddefn{#1}}
\enddefn{\relax}

\newenvironment{Defn}[1][]{\refstepcounter{defn}%
\list{\relax}{\topsep=5pt \settowidth{\labelwidth}{0pt}
\leftmargin=0pt \labelsep=0pt \rightmargin\leftmargin\itemsep=0pt
\parsep=0pt\itemindent=0pt\listparindent=\parindent}
\item[] {\scshape Definition\ \thedefn\@periodafterdefn\
\ifx#1\undefiend\relax\else\@llbrac #1\@rrbrac\@enddefn\fi}\par\vskip0pt\noindent\ignorespaces}
{\endlist}

\newenvironment{fldefn}[1][]{\refstepcounter{defn}\list{\relax}{\topsep=5pt
\settowidth{\labelwidth}{0pt} \leftmargin=0pt
\labelsep=0pt\rightmargin\leftmargin\itemsep=0pt
\parsep=0pt\itemindent=0pt\listparindent=\parindent}
\item[] {\scshape Definition\ \thedefn\@periodafterdefn\ \ifx#1\undefiend
\else\@llbrac #1\@rrbrac\@enddefn\fi}\hskip.3em} {\endlist}

\newenvironment{flthem}[1][]{\refstepcounter{defn}\list{\relax}{\topsep=3pt
\settowidth{\labelwidth}{0pt} \leftmargin=0pt
\thmfont\labelsep=0pt\rightmargin\leftmargin\itemsep=0pt
\parsep=0pt\itemindent=0pt\listparindent=\parindent}
\item[] {\bf Theorem\ \thedefn\@periodafterdefn\ \ifx#1\undefiend
\else\@llbrac #1\@rrbrac\@enddefn\fi}\hskip.3em\itshape} {\endlist}

\newcounter{propn}
\def\thepropn{\arabic{propn}}
\def\periodafterpropn#1{\gdef\@periodafterpropn{#1}}
\periodafterpropn{.}
\def\periodafterlem#1{\gdef\@periodafterlem{#1}}
\periodafterlem{.}

\def\llbrac#1{\gdef\@llbrac{#1}}
\def\rrbrac#1{\gdef\@rrbrac{#1}}
\llbrac{\text{\rm (}}
\rrbrac{\text{\rm )}}

\newenvironment{Prop}[1][]{\refstepcounter{propn}%
\list{\relax}{\topsep=5pt \settowidth{\labelwidth}{0pt}
\leftmargin=0pt \labelsep=0pt \rightmargin\leftmargin\itemsep=0pt
\parsep=0pt\itemindent=0pt\listparindent=\parindent}
\item[] {\scshape Proposition\ \thepropn\@periodafterpropn\
\ifx#1\undefiend\relax\else\@llbrac #1\@rrbrac\fi}\par\vskip0pt\noindent\ignorespaces}
{\endlist}

\newenvironment{flProp}[1][]{\refstepcounter{propn}\list{\relax}{\topsep=5pt
\settowidth{\labelwidth}{0pt} \leftmargin=0pt
\labelsep=0pt\rightmargin\leftmargin\itemsep=0pt
\parsep=0pt\itemindent=0pt\listparindent=\parindent}
\item[] {\scshape \bf Proposition\ \thepropn\@periodafterpropn\ \ifx#1\undefiend
\else\@llbrac #1\@rrbrac\fi}\hskip.3em\itshape} {\endlist}

\newcounter{Lem}
\def\theLem{\arabic{Lem}}

\newenvironment{Lem}[1][]{\refstepcounter{Lem}%
\list{\relax}{\topsep=5pt \settowidth{\labelwidth}{0pt}
\leftmargin=0pt \labelsep=0pt \rightmargin\leftmargin\itemsep=0pt
\parsep=0pt\itemindent=0pt\listparindent=\parindent}
\item[] {\scshape Lemma\ \theLem\@periodafterlem\
\ifx#1\undefiend\relax\else\@llbrac #1\@rrbrac\fi}\par\vskip0pt\noindent\ignorespaces}
{\endlist}

%\newenvironment{flLem}[1][]{\refstepcounter{Lem}\list{\relax}{\topsep=5pt
%\settowidth{\labelwidth}{0pt} \leftmargin=0pt
%\labelsep=0pt\rightmargin\leftmargin\itemsep=0pt
%\parsep=0pt\itemindent=0pt\listparindent=\parindent}
%\item[] {\scshape\bf Lemma\ \thepropn\@periodafterlem\ \ifx#1\undefiend
%\else\@llbrac #1\@rrbrac\fi}\hskip.3em} {\endlist}

\newenvironment{flLem}[1][]{\refstepcounter{propn}\list{\relax}{\topsep=5pt
\settowidth{\labelwidth}{0pt} \leftmargin=0pt
\labelsep=0pt\rightmargin\leftmargin\itemsep=0pt
\parsep=0pt\itemindent=0pt\listparindent=\parindent}
\item[] {\scshape\bf Lemma\ \thepropn\@periodafterlem\ \ifx#1\undefiend
\else\@llbrac #1\@rrbrac\fi}\hskip.3em\itshape} {\endlist}


\newenvironment{flPropl}[1][]{\refstepcounter{propn}\list{\relax}{\topsep=5pt
\settowidth{\labelwidth}{0pt} \leftmargin=0pt
\labelsep=0pt\rightmargin\leftmargin\itemsep=0pt
\parsep=0pt\itemindent=0pt\listparindent=\parindent}
\item[] {\scshape Problem\ \thepropn\@periodafterpropn\ \ifx#1\undefiend
\else\@llbrac #1\@rrbrac\fi}\hskip.3em\itshape} {\endlist}


\newenvironment{flcase}[1][]{\refstepcounter{propn}\list{\relax}{\topsep=5pt
\settowidth{\labelwidth}{0pt} \leftmargin=0pt
\labelsep=0pt\rightmargin\leftmargin\itemsep=0pt
\parsep=0pt\itemindent=0pt\listparindent=\parindent}
\item[] {\scshape Case\ \thepropn\@periodafterpropn\ \ifx#1\undefiend
\else\@llbrac #1\@rrbrac\fi}\hskip.3em} {\endlist}

\newcounter{rem}
\def\theLem{\arabic{rem}}

\newenvironment{flrem}[1][]{\refstepcounter{propn}\list{\relax}{\topsep=5pt
\settowidth{\labelwidth}{0pt} \leftmargin=0pt
\labelsep=0pt\rightmargin\leftmargin\itemsep=0pt
\parsep=0pt\itemindent=0pt\listparindent=\parindent}
\item[] {\scshape Remark\ \thepropn\@periodafterpropn\ \ifx#1\undefiend
\else\@llbrac #1\@rrbrac\fi}\hskip.3em} {\endlist}

%\newcounter{defin}
%\def\thedefn{\arabic{defin}}
%\def\periodafterdefn#1{\gdef\@periodafterdefin{#1}}
%\periodafterdefn{.}
%
%\newenvironment{Defin}[1][]{\refstepcounter{defin}\list{\relax}{\topsep=5pt
%\settowidth{\labelwidth}{0pt} \leftmargin=0pt
%\labelsep=0pt\rightmargin\leftmargin\itemsep=0pt
%\parsep=0pt\itemindent=0pt\listparindent=\parindent}
%\item[] {\scshape Definition\ \thedefn\@periodafterdefn\ \ifx#1\undefiend
%\else(#1)\fi}\vskip0pt} {\endlist}

\def\proofname{\textit{Proof}}

\renewenvironment{proof}[1][\proofname]{\par
  \normalfont \topsep6\p@\@plus6\p@\relax
  \labelsep 0.5em
  \trivlist
  \item[\hskip\labelsep\hskip1em\textsc{#1}.~~]\ignorespaces
}{{\hfill$\square$}\endtrivlist\@endpefalse}

\def\editorialhead#1{{\fontfamily{\rmdefault}\fontsize{21}{21}\fontseries{b}\fontshape{it}\selectfont\leftskip0pt\rightskip0pt
plus1fill \noindent #1\par\vskip\baselineskip}}

\makeindex

\newcommand\abc[1]{\par{\addvspace{6pt plus 2pt minus 2pt}}\noindent\textit{#1}:\ }
\newcommand\abcd[1]{\par{\addvspace{6pt plus 2pt minus 2pt}}\noindent\textrm{#1}:\ }

%\newcommand{\mcc}[1]{\multicolumn{1}{c}{#1}}
%\newcommand{\mcl}[1]{\multicolumn{1}{l}{#1}}
%\newcommand{\mcr}[1]{\multicolumn{1}{r}{#1}}

\def\au#1{\index{#1}}

\def\mcc#1{\multicolumn{1}{c}{#1}}
%%%%%%%%%%Funding and Conflict%%%

\def\textcolon{\text{\rm :}}

 \def\pparagraph{%
   \@startsection{paragraph}{4}{\z@}
   {-6\p@}
   {-.4em}
   {\reset@font\itshape}}

%\font\txfont=CRONOS at 10pt
%\def\ndash{\mbox{\txfont\char"96}}
%\def\mdash{\mbox{\txfont\char"97}}
%
%\font\sltxfont=CRONOS at 6pt
%\def\slndash{\mbox{\sltxfont\char"96}}
%\def\slmdash{\mbox{\sltxfont\char"97}}

\endinput

\interlinepenalty -100000
\penalty -30000
\newcount\linenopenalty\linenopenalty=-100000
\def\linenopenalty{-100000\relax}
\mathchardef\linenopenaltypar=32000

\def\@outputpage{%
\begingroup           % the \endgroup is put in by \aftergroup
  \let \protect \noexpand
  \@resetactivechars
  \@parboxrestore
  \shipout \vbox{%
    \set@typeset@protect
    \aftergroup \endgroup
    \aftergroup \set@typeset@protect
                                % correct? or just restore by ending
                                % the group?
  \if@specialpage
    \global\@specialpagefalse\@nameuse{ps@\@specialstyle}%
  \fi
  \if@twoside
    \ifodd\count\z@ \let\@thehead\@oddhead \let\@thefoot\@oddfoot
         \let\@themargin\oddsidemargin
    \else \let\@thehead\@evenhead
       \let\@thefoot\@evenfoot \let\@themargin\evensidemargin
    \fi
  \fi
  \reset@font
  \normalsize
  \normalsfcodes
  \let\label\@gobble
  \let\index\@gobble
  \let\glossary\@gobble
  \baselineskip\z@skip \lineskip\z@skip \lineskiplimit\z@
    \@begindvi
    \vskip \topmargin
    \moveright\@themargin \vbox {%
      \setbox\@tempboxa \vbox to\headheight{%
        \vfil
        \color@hbox
          \normalcolor
          {\hb@xt@\textwidth{\@thehead}}%
        \color@endbox
        }%                        %% 22 Feb 87
      \dp\@tempboxa \z@
      \box\@tempboxa
      \vskip \headsep
      {\box\@outputbox}
      \baselineskip \footskip
      \color@hbox
        \normalcolor
        {\hb@xt@\textwidth{\@thefoot}}%
      \color@endbox
      }%
    }%
  \global\let\if@newlist\@@if@newlist
  \global \@colht \textheight
  \stepcounter{page}%
  %\ifnum\count0=341\thispagestyle{empty}\fi
  %\ifnum\count0=61\thispagestyle{empty}\fi
  \let\firstmark\botmark
}

% \def\sparagraph{%
%   \@startsection{paragraph}{4}{\z@}
%   {-6\p@}
%   {-.4em}
%   {\reset@font\itshape}}

%%
%% End of file `bioinfo.cls'.

