from stytra import Stytra
from stytra.stimulation.stimuli import (
    MovingGratingStimulus,
    GratingStimulus,
    PositionStimulus,
    BackgroundStimulus,
)

from PyQt5.QtGui import QTransform

from stytra.stimulation.stimuli.conditional import TwoRadiusCenteringWrapper
from stytra.stimulation import Protocol
from lightparam import Param
import pandas as pd
import numpy as np
from stytra_config import ConfiguredStytra


class FishRelativeRotationOnlyStimulus(BackgroundStimulus):
    def __init__(self, *args, theta_change_threshold=12 * np.pi / 180, **kwargs):
        super().__init__(*args, **kwargs)
        self.stored_fish_theta = None
        self.fish_theta_change_threshold = theta_change_threshold

    @property
    def theta_total(self):
        if self.stored_fish_theta is not None:
            return self.theta + self.stored_fish_theta - np.pi / 2
        return self.theta

    def get_transform(self, w, h, x, y):
        _, _, theta_fish = self._experiment.estimator.get_position() # underscore means that this will be discarded
        if self.stored_fish_theta is None or \
                np.abs(np.mod(theta_fish-self.stored_fish_theta, 2*np.pi)) > self.fish_theta_change_threshold:
            self.stored_fish_theta = theta_fish

        rot_fish = (self.stored_fish_theta - np.pi / 2) * 180 / np.pi
        xc = w/2
        yc = h/2
        return super().get_transform(w, h, x, y) * (QTransform().translate(xc, yc).rotate(rot_fish).translate(-xc, -yc))


class GratingsTrackingStimulus(FishRelativeRotationOnlyStimulus, MovingGratingStimulus):
    pass


class ColorOmrFreelySwimmingProtocol(Protocol):
    name = "free_omr_lr_cl_contrasts_v02"

    stytra_config = dict(
        display=dict(min_framerate=50),
        tracking=dict(method="fish", embedded=False, estimator="position"),

    )

    def __init__(self):
        super().__init__()

        self.inter_stim_pause = Param(10., limits=(0, 300))
        self.grating_vel = Param(10., limits=(-50, 50))
        self.grating_duration = Param(15., limits=(0, 300))
        self.grating_cycle = Param(10, limits=(0, 300))  # spatial period of the grating
        self.n_trials = Param(5, (0, 100))
        #self.contrast = Param(10, (0, 150))
        #self.n_rep_internal = Param(10, limits=(0, 300))  # nr of total red and green alternations (10 each; 20 total)


    def get_stim_sequence(self):
            stimuli = []

            v = self.grating_vel
            d = self.grating_duration
            p = self.inter_stim_pause

            t = [0, p, p,  p+d, p+d, 2*p+d, 2*p+d, 2*p+2*d]
            vel = [0, 0, v, v, 0, 0, -v, -v]
            contrast = [55, 60, 65, 70, 90]

            df = pd.DataFrame(dict(t=t, vel_x=vel))
            for i in range(self.n_trials):
                stimuli.append(TwoRadiusCenteringWrapper(stimulus=GratingsTrackingStimulus(
                        df_param=df,
                        grating_period=self.grating_cycle,
                        grating_col_1=(contrast[i],)*3,
                        #grating_col_1= (self.contrast,)*3,
                        grating_angle=0,
                    ),  r_out=45, r_in=40
                ))

            return stimuli



if __name__ == "__main__":
    s = Stytra(protocol=ColorOmrFreelySwimmingProtocol())