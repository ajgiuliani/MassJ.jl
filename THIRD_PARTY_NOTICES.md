# Third-party notices

MassJ.jl as a whole is distributed under the GNU General Public License
v3.0 or later (see [LICENSE](LICENSE)). It incorporates ideas and techniques
from the following third-party works, whose original licenses and notices are
reproduced below as required.

---

## julia_mzML_imzML

The fast DOM-free streaming imzML reader in `src/imzml_stream.jl` adapts the
parsing strategy of **julia_mzML_imzML** by Dr Robert Winkler — reading the
`<spectrumList>` as a text stream rather than building a full XML DOM, and
reading binary arrays directly from the companion `.ibd` file.

- Project: <https://codeberg.org/LabABI/julia_mzML_imzML>
- Associated publication:
  I. Rosas-Román, H. Guillén-Alonso, A. Moreno-Pedraza, R. Winkler,
  *Anal. Chem.* **2024**. DOI:
  [10.1021/acs.analchem.3c05853](https://doi.org/10.1021/acs.analchem.3c05853)
- Test/benchmark data: R. Winkler (2023), Zenodo,
  DOI: [10.5281/zenodo.10084132](https://doi.org/10.5281/zenodo.10084132)

License (MIT):

```
MIT License

Copyright (c) 2023 Dr Robert Winkler

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
