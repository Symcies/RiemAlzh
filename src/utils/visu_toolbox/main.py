import argparse

from utils import ui_management

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
    ui_management.plot(args['p'], args['i'], args['o'], args['t'])

if __name__ == '__main__':
    main()