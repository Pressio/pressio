
import sys
import re

new_header = '''/*
//@HEADER
// ************************************************************************
//
// {:^72}
//                         dharma_mockup
//              Copyright (C) 2015 Sandia Corporation
// This file was adapted from its original form in the tinympl library.
// The original file bore the following copyright:
//   Copyright (C) 2013, Ennio Barbaro.
// See LEGAL.md for more information.
//
// Under the terms of Contract DE-NA-0003525 with NTESS, LLC,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact somebody@sandia.gov
//
// ************************************************************************
//@HEADER
*/

'''

for fname in sys.argv[1:]:
  with open(fname) as f:
    contents = f.read()
    lines = contents.splitlines()
    if re.search(
      r'Use, modification, and distribution is subject to the Boost Software',
      "\n".join(lines[:11])
    ):
      rest_of_file = "\n".join(lines[11:])
    else:
      rest_of_file = None
  if rest_of_file is not None:
    with open(fname, "w+") as f:
      f.write(new_header.format(fname.split("/")[-1]))
      f.write(rest_of_file)


for fname in sys.argv[1:]:
  with open(fname) as f:
    contents = f.read()
    if 'variadic/' not in fname:
      contents = re.sub(
        r'#include <tinympl/(\S+)>',
        r'#include "\1"',
        contents
      )
    else:
      contents = re.sub(
        r'#include <tinympl/variadic/(\S+)>',
        r'#include "\1"',
        contents
      )
      contents = re.sub(
        r'#include <tinympl/(\S+)>',
        r'#include "../\1"',
        contents
      )
      content_lines = contents.splitlines()
      if content_lines[0].strip() == '':
        content_lines = content_lines[1:]
      contents = "\n".join(content_lines)
  with open(fname, "w+") as f:
    f.write(contents)


