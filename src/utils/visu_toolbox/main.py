import argparse

from visu import Visu

def main():
    parser = argparse.ArgumentParser(description='Displays the curves for the different individuals')
    parser.add_argument('-p', action='store',
                        help='path of the final population file')
    parser.add_argument('-i', action='store',
                        help='path of the final individual file')
    parser.add_argument('-o', action='store',
                        help='path of the final observations file')
    parser.add_argument('-t', action='store',
                        help='type of the model')

    args = vars(parser.parse_args())
    visu = Visu(args['t'], args['p'], args['i'], args['o'])
    visu.init_model()
    visu.plot_mean()
    visu.plot_patients([10, 20], True)
    visu.model.plot_stand_dev("tau")
    visu.model.plot_stand_dev("ksi")
    visu.hold_plot()


if __name__ == '__main__':
    main()