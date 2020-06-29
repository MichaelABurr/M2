###############################################################################
## This directory contains the rules for generating grammar files for editors
## and syntax highlighing engines as well as the M2-mode package for Emacs.
## - M2-editors generates grammar files
## - M2-emacs   generates the M2-mode package for Emacs

# TODO: generate vim grammar
# TODO: generate textmate grammar

set(M2 ${M2_DIST_PREFIX}/${M2_INSTALL_BINDIR}/M2)
set(M2_ARGS --script)

###############################################################################
## Generating the grammar files is mainly done in make-M2-symbols.m2.
## To add to the M2-editors target, add grammar filenames here and in .gitignore.
set(GRAMMAR_FILES
  atom/macaulay2.cson
  prism/macaulay2.js
#  rouge/macaulay2.rb
  vim/m2.vim.syntax
  vim/m2.vim.dict
  emacs/M2-symbols.el)

add_custom_target(M2-editors DEPENDS ${GRAMMAR_FILES})
add_custom_command(OUTPUT ${GRAMMAR_FILES}
  COMMENT "Generating M2 grammar files for editors"
  COMMAND ${M2} ${M2_ARGS} ${CMAKE_CURRENT_SOURCE_DIR}/make-M2-symbols.m2
  MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/make-M2-symbols.m2
  DEPENDS M2-core)

###############################################################################
## Generating the M2-mode package for Emacs

# TODO: only the emacs target is installed in usr-dist by default
set(M2_INSTALL_EMACSDIR
  ${M2_DIST_PREFIX}/${M2_INSTALL_DATAROOTDIR}/emacs/site-lisp/Macaulay2)
file(COPY emacs/ DESTINATION ${M2_INSTALL_EMACSDIR}
  FILES_MATCHING PATTERN "M2*" PATTERN "*.in" EXCLUDE)

add_custom_target(M2-emacs ALL
  COMMENT "Generating Emacs package"
  COMMAND ${M2} ${M2_ARGS} ${CMAKE_CURRENT_SOURCE_DIR}/emacs/make-M2-emacs-help.m2
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/emacs/M2-symbols.el ${M2_INSTALL_EMACSDIR}
  COMMAND gzip -nf9 M2-symbols.el
  COMMAND ${CMAKE_COMMAND} -E copy M2-emacs-help.txt M2-emacs.m2 M2-symbols.el.gz ${CMAKE_CURRENT_SOURCE_DIR}/emacs
  WORKING_DIRECTORY ${M2_INSTALL_EMACSDIR}
  DEPENDS M2-editors)