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

from chimerax.mouse_modes import MouseMode


class DrawBinderMouseMode(MouseMode):
    name = 'draw binder'
    icon_file = 'draw_binder.png'

    def __init__(self, session):
        MouseMode.__init__(self, session)
        self._bond_rot = None
        self._marker_num = 0
        self._marker_coords = []
        self._last_x = None
        self._last_y = None
        self._last_coord = None
        self._speed_factor = 2
        self.element_counts = {}
        # Degrees.  Only applies to drag with 3d pointer.
        self._minimum_angle_step = 2

    def mouse_down(self, event):
        MouseMode.mouse_down(self, event)
        from chimerax.core.commands import run
        # for i in range(1, 1+self._marker_num):
        from chimerax.markers import MarkerSet
        mlist = [m for m in self.session.models.list(
            type=MarkerSet) if m.id == (2,)]
        if mlist:
            self._marker_num = max(mlist[0].residues.numbers)+1
        else:
            self._marker_num = 1
        # run(self.session, 'marker delete #2')

        x, y = event.position()
        coord = self.window_scene_position(x, y)
        # self._marker_coords = [coord]
        if coord is None:
            if mlist:
                self._marker_num = max(mlist[0].residues.numbers)
            else:
                self._marker_num = 0
            # self._marker_coords = []
            return

        # marker #2 position 102.6,137.5,78.41 color white radius 1
        run(self.session, 'marker #2 position %f,%f,%f color white radius 0.2' %
            (coord[0], coord[1], coord[2]))
        if self._marker_num > 1:
            run(self.session, 'marker link #2:%d,%d color white radius 0.2' %
                (self._marker_num-1, self._marker_num))
        self._last_coord = coord
        self._last_x, self._last_y = event.position()
        # pick = self._picked_bond(event)
        # self._bond_rot = self._picked_bond_rotation(
        #     pick, move_smaller_side=not event.shift_down())
        # self._first_coord = self.window_scene_position(event)

    def able_draw(self, x, y, threshold=1):
        # coord = self.window_scene_position(x, y)
        if self._last_x is None:
            return False
        dx = x - self._last_x
        dy = y - self._last_y
        if dx*dx + dy*dy > threshold * threshold:

            return True
        else:
            return False

    def spline(self, x, y, threshold=3, base=2.5):
        from math import sqrt
        from chimerax.core.commands import run
        coord = self.window_scene_position(x, y)
        if coord is None:
            return
        distance = sqrt((coord[0]-self._last_coord[0])**2 +
                        (coord[1]-self._last_coord[1])**2 + (coord[2]-self._last_coord[2])**2)
        if distance > threshold:
            # spline_count = int(distance//base)
            spline_count = 1
            self.session.logger.status(" spline_count: %d" % (spline_count))
        else:
            spline_count = 0

        for i in range(1, spline_count+1):
            # print("spline_count: %d" % (spline_count))
            spline_x = x-(x-self._last_x)*i/(spline_count+1)
            spline_y = y-(y-self._last_y)*i/(spline_count+1)
            spline_coord = self.window_scene_position(spline_x, spline_y)
            if spline_coord is None:
                continue
            self._marker_num += 1
            # self._marker_coords.append(spline_coord)
            run(self.session, 'marker #2 position %f,%f,%f color white radius 0.2' % (
                spline_coord[0],
                spline_coord[1],
                spline_coord[2]))

            # marker link #2:10,2 color #659cef radius 0.2
            if self._marker_num > 1:
                run(self.session, 'marker link #2:%d,%d color white radius 0.2' %
                    (self._marker_num-1, self._marker_num))
        self._marker_num += 1
        # self._marker_coords.append(coord)
        run(self.session, 'marker #2 position %f,%f,%f color white radius 0.2' %
            (coord[0], coord[1], coord[2]))
        if self._marker_num > 1:
            run(self.session, 'marker link #2:%d,%d color white radius 0.2' %
                (self._marker_num-1, self._marker_num))
        self._last_coord = coord
        self._last_x = x
        self._last_y = y

    def mouse_drag(self, event):
        dx, dy = self.mouse_motion(event)
        x, y = event.position()
        if self.able_draw(x, y, 25):

            self.spline(x, y, threshold=2.5)

            # from chimerax.core.commands import run
            # run(self.session, 'build start atom #1 pos %f,%f,%f' %
            #     (coord[0], coord[1], coord[2]))
            # run(self.session, 'bond /het:%d@He /het:%d@He reasonable false' %
            #     (self._marker_num-1, self._marker_num))
            # self.session.logger.status(" distance: %f" % (distance))
    def save_xyz(self, path):
        """Write an XYZ file from given models, or all models if None.
        """
        # Open path with proper encoding; 'open_output' automatically
        # handles compression if the file name also has a compression
        # suffix (e.g. .gz)
        from chimerax.io import open_output
        # f = open_output(path, 'utf-8')
        import numpy as np
        # np.save(path, self._marker_coords)
        # np.savetxt(
        #     path,
        #     self._marker_coords,
        #     fmt="%.3f",
        #     delimiter=",",
        # )
        # for i in range(len(self._marker_coords)):
        #     c = self._marker_coords[i]
        #     print("%.3f %.3f %.3f" %
        #           (c[0], c[1], c[2]), file=f)
        # f.close()

        # Notify user that file was saved
        self.session.logger.info("Saved binder containing %d atoms to %s"
                                 % (self._marker_num, path))

    def window_scene_position(self, win_x, win_y):

        ses = self.session
        view = ses.main_view
        xyz_min, xyz_max = view.clip_plane_points(win_x, win_y)
        p = ses.main_view.picked_object_on_segment(
            xyz_min, xyz_max, max_transparent_layers=0)
        if p is None:
            return None
        else:
            scene_pos = .05*xyz_min + .95*p.position
            # from chimerax.core.commands import log_equivalent_command
            # distance = scene_pos - xyz_min
            # log_equivalent_command(self.session, 'build start atom #1 pos %f,%f,%f' % (
            #     coord[0], coord[1], coord[2]))
        # scene_point = (0,0,0) if xyz_min is None else (.9*xyz_min + .1*xyz_max)

        # # Convert camera coordinates to scene coordinates
        # rot = view.camera.position.zero_translation()
        # from chimerax.geometry import translation
        # scene_pos = translation(scene_point) * rot
        # coord = scene_pos.origin()
        return scene_pos

    def mouse_up(self, event):
        from chimerax.core.commands import log_equivalent_command
        MouseMode.mouse_up(self, event)
        import os
        pwd = os.getcwd()
        save_file = os.path.join(pwd, "binder.npy")
        if save_file is not None:
            self.save_xyz(save_file)
        # coord = self.window_scene_position(event)
        # for i in range(1, self._marker_num):
        #     from chimerax.core.commands import run
        #     run(self.session, 'bond /het:%d@He /het:%d@He reasonable false' % (i, i+1))
        # log_equivalent_command(self.session, 'draw %d %d %d' %
        #                        (coord[0], coord[1], coord[2]))
        # initial_coord = self._first_coord
        # log_equivalent_command(self.session, 'draw %d %d %d' % (
        #     initial_coord[0], initial_coord[1], initial_coord[2]))
        # from chimerax.atomic import AtomicStructure
        # if coord is not None:
        #     for m in self.session.models:
        #         if isinstance(m, AtomicStructure):
        #             # log_equivalent_command(self.session, '%s'% m)
        #             # log_equivalent_command(self.session, 'build start atom %s pos %f,%f,%f' % (
        #             #     m, coord[0], coord[1], coord[2]))
        #             # log_equivalent_command(self.session, 'build start atom %s pos %f,%f,%f' % (
        #             #     m, initial_coord[0], initial_coord[1], initial_coord[2]))
        #             from chimerax.core.commands import run
        #             run(self.session, 'build start atom #1 pos %f,%f,%f' %
        #                 (coord[0], coord[1], coord[2]))
        # run(self.session, 'build start atom #1 pos %f,%f,%f' %
        #     (coord_max[0], coord_max[1], coord_max[2]))
        # log_equivalent_command(self.session, 'build start atom #1 pos %f,%f,%f' % (
        #     coord[0], coord[1], coord[2]))
        # scene = self.scene()
        # self._drag_box = scene.addRect(1, 2, 3, 4)
        # self._log_command()
        # self._delete_bond_rotation()

    def _log_command(self):
        br = self._bond_rotation
        if br:
            log_torsion_command(br)

    def wheel(self, event):
        pick = self._picked_bond(event)
        br = self._picked_bond_rotation(
            pick, move_smaller_side=not event.shift_down())
        if br:
            d = event.wheel_value()
            br.angle += d
            self.session.bond_rotations.delete_rotation(br)

    def _picked_bond(self, event):
        x, y = event.position()
        pick = self.session.main_view.picked_object(x, y)
        return pick

    def _picked_bond_rotation(self, pick, move_smaller_side=True):
        from chimerax.atomic import PickedBond
        if isinstance(pick, PickedBond):
            # from .manager import BondRotationError
            # try:
            #     br = self.session.bond_rotations.new_rotation(
            #         pick.bond, move_smaller_side=move_smaller_side)
            #     self.session.logger.status('Rotating bond %s' % str(pick.bond))
            # except BondRotationError as e:
            #     self.session.logger.status(str(e))
            self.session.logger.status('Rotating bond %s' % str(pick.bond))
            br = None
        else:
            br = None
        return br

    @ property
    def _bond_rotation(self):
        br = self._bond_rot
        if br and br.bond.deleted:
            self._bond_rot = br = None
        return br

    def _delete_bond_rotation(self):
        br = self._bond_rotation
        if br is not None:
            self.session.bond_rotations.delete_rotation(br)
            self._bond_rot = None

    def vr_press(self, event):
        # Virtual reality hand controller button press.
        pick = event.picked_object(self.view)
        self._bond_rot = br = self._picked_bond_rotation(pick)
        if br:
            br.bond.selected = True

        # Move the side of the bond the VR click is closest to.
        # Would like to have a command to enable this mode for rotating bonds
        # with small ligands
        move_closer_side = False
        if move_closer_side and br is not None:
            atom1 = br.moving_side
            atom2 = br.bond.other_atom(atom1)
            p = event.tip_position
            from chimerax.geometry import distance
            if distance(p, atom2.scene_coord) < distance(p, atom1.scene_coord):
                br.moving_side = atom2

    def vr_motion(self, event):
        # Virtual reality hand controller motion.
        br = self._bond_rotation
        if br:
            axis, angle = event.motion.rotation_axis_and_angle()
            from chimerax.geometry import inner_product
            if inner_product(axis, br.axis) < 0:
                angle = -angle
            angle_change = self._speed_factor * angle
            if abs(angle_change) < self._minimum_angle_step:
                return "accumulate drag"
            br.angle += angle_change

    def vr_release(self, event):
        # Virtual reality hand controller button release.
        br = self._bond_rotation
        if br:
            br.bond.selected = False
            self._log_command()
            self._delete_bond_rotation()


def log_torsion_command(bond_rotator):
    bond = bond_rotator.rotation.bond
    ms_atom = bond_rotator.moving_side
    fs_atom = bond.other_atom(ms_atom)
    ms_atom2 = _connected_atom(ms_atom, fs_atom)
    fs_atom2 = _connected_atom(fs_atom, ms_atom)
    if ms_atom2 is None or fs_atom2 is None:
        return 		# No connected atom to define a torsion
    side = '' if bond.smaller_side is ms_atom else 'move large'
    from chimerax.geometry import dihedral
    torsion = dihedral(fs_atom2.scene_coord, fs_atom.scene_coord,
                       ms_atom.scene_coord, ms_atom2.scene_coord)

    atom_specs = '%s %s %s %s' % (fs_atom2.string(style='command'), fs_atom.string(style='command'),
                                  ms_atom.string(style='command'), ms_atom2.string(style='command'))

    # Use simpler atom spec for the common case of rotating a side chain.
    res = ms_atom.residue
    if ms_atom2.residue is res and fs_atom.residue is res and fs_atom2.residue is res:
        # serial_number indicates a duplicate atom name.
        if 'serial_number' not in atom_specs:
            atom_specs = '%s@%s,%s,%s,%s' % (res.string(style='command'),
                                             ms_atom2.name, ms_atom.name, fs_atom.name, fs_atom2.name)

    cmd = 'torsion %s %.2f %s' % (atom_specs, torsion, side)
    ses = ms_atom.structure.session
    from chimerax.core.commands import run
    run(ses, cmd)


def _connected_atom(atom, exclude_atom):
    for a in atom.neighbors:
        if a is not exclude_atom:
            return a
    return None
