'''
Author: error: error: git config user.name & please set dead value or install git && error: git config user.email & please set dead value or install git & please set dead value or install git
Date: 2023-12-19 19:40:18
LastEditors: error: error: git config user.name & please set dead value or install git && error: git config user.email & please set dead value or install git & please set dead value or install git
LastEditTime: 2024-07-04 15:52:37
FilePath: \chimerax\draw_line\src\__init__.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''
# vim: set expandtab shiftwidth=4 softtabstop=4:

# === UCSF ChimeraX Copyright ===
# Copyright 2022 Regents of the University of California. All rights reserved.
# The ChimeraX application is provided pursuant to the ChimeraX license
# agreement, which covers academic and commercial uses. For more details, see
# <http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html>
#
# This particular file is part of the ChimeraX library. You can also
# redistribute and/or modify it under the terms of the GNU Lesser General
# Public License version 2.1 as published by the Free Software Foundation.
# For more details, see
# <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>
#
# THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
# EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. ADDITIONAL LIABILITY
# LIMITATIONS ARE DESCRIBED IN THE GNU LESSER GENERAL PUBLIC LICENSE
# VERSION 2.1
#
# This notice must be embedded in or attached to all copies, including partial
# copies, of the software or any revisions or derivations thereof.
# === UCSF ChimeraX Copyright ===

from chimerax.core.toolshed import BundleAPI


class _SaveCurveAPI(BundleAPI):

    @staticmethod
    def initialize(session, bundle_info):
        """Install bond-rotation manager into existing session"""
        # from .manager import BondRotationManager
        # session.draw_binders = BondRotationManager(session, bundle_info)
        if session.ui.is_gui:
            mm = session.ui.mouse_modes
            from .save_curve import SaveCurveMouseMode
            mm.add_mode(SaveCurveMouseMode(session))

    # @staticmethod
    # def finish(session, bundle_info):
    #     """De-install bond-rotation manager from existing session"""
    #     del session.draw_binders

    @staticmethod
    def register_command(command_name, logger):
        # 'register_command' is lazily called when the command is referenced
        from . import cmd
        cmd.register_command(command_name, logger)


bundle_api = _SaveCurveAPI()