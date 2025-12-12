import numpy as np
import matplotlib.pyplot as plt

# Free Streaming Model
# particles: Human-defined plasma with e and i (m = 100 )


class Setting:
    def __init__(
        self, electronTime=True, simpleFSM=True, ThermalSpeed=True, ParticleFlux=True
    ):
        self.sigma = 10
        self.electronTime = electronTime
        self.simpleFSM = simpleFSM
        self.ThermalSpeed = ThermalSpeed
        self.ParticleFlux = ParticleFlux


class Ion:
    def __init__(self, name, mass, charge):
        self.name = name
        self.mass = mass
        self.charge = charge


class FSM:
    def __init__(self, setting, name, ions):
        self.name = name
        self.ions = ions
        self.setting = setting

    def introduce(self):
        print(
            f" I'm {self.name} plasma with ions {self.ions.name} (mass = {self.ions.mass} and charge = {self.ions.charge})"
        )

    def Output(self, time):
        t = time  # transition time denfined by sigma/sound speed
        sigma = self.setting.sigma

        if self.setting.simpleFSM:
            if self.setting.ParticleFlux:
                flux = sigma / (t**2) * np.exp(-0.5 * (sigma**2) / (t**2))
            else:
                flux = (
                    sigma
                    / (t**2)
                    * np.exp(-0.5 * (sigma**2) / (t**2))
                    * (sigma**2 / (t**2) + 1)
                )
        else:
            if self.setting.ParticleFlux:
                flux = (
                    sigma
                    * t
                    / (1 + t**2) ** (3 / 2)
                    * np.exp(-0.5 * sigma**2 / (1 + t**2))
                )
            else:
                flux = (
                    sigma
                    * t
                    / (1 + t**2) ** (3 / 2)
                    * np.exp(-0.5 * sigma**2 / (1 + t**2))
                ) * (sigma**2 * t**2 / (1 + t**2) ** 2 + 1.5 / (1 + t**2) + 1)

        if self.setting.ThermalSpeed:
            flux = flux * np.sqrt(2 / self.ions.mass)
        else:
            flux = flux

        return flux

    def TimeTransfer(self, time):

        if self.setting.electronTime:
            t = (
                time * np.sqrt(0.5 * self.ions.mass) / self.setting.sigma
            )  # transfer into e- transition time
        else:
            t = (
                time
                * np.sqrt(0.5 * self.ions.mass)
                / self.setting.sigma
                / np.sqrt(1836)
            )  # H+ transition time
        return t

    def show(self, time):

        f = self.Output(time)
        t = self.TimeTransfer(time)

        plt.figure(figsize=(8, 4))

        if self.setting.electronTime:
            plt.xlabel(r"Time/$\tau_e$")
        else:
            plt.xlabel(r"Time/$\tau_H^+$")

        if self.setting.ThermalSpeed:
            plt.ylabel(r"Flux($n_ev^{th}_e$)")
        else:
            plt.ylabel(r"Flux($n_ec_s$)")

        plt.plot(t, f, label="f(t)", color="blue")
        # plt.xlim(0, 4)
        plt.title("Electron Transition Flux vs Time")
        plt.grid(True)
        plt.show()


if __name__ == "__main__":
    fsm = FSM(
        Setting(True, False, True, True), "Human-defined Model", Ion("H+", 1836, 1)
    )

    time = np.linspace(0.01, 40, 1000)

    fsm.show(time)
