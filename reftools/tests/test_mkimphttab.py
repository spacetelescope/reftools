import os
import tempfile

from .. import mkimphttab


# make sure that when making an imphttab there are no errors
def test_mkimphttab():
  tempdir = tempfile.gettempdir()

  output = os.path.join(tempdir, 'test_out_imp.fits')

  base_mode = 'acs,sbc'

  useafter='DUMMY'

  detector='sbc'

  try:
    mkimphttab.create_table(output, base_mode, detector, useafter)

  except:
    if os.path.exists(output):
      os.remove(output)

    raise

  else:
    assert os.path.exists(output)
    os.remove(output)
