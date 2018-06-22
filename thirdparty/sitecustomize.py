import site

import os
LOCAL_PATH = os.path.dirname(os.path.realpath(__file__))

SITE_PACKAGES_PATH = os.path.join(LOCAL_PATH, 'site-packages')

site.addsitedir(SITE_PACKAGES_PATH)
