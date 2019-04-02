import json
import random
import sys
import getopt
import uuid
import datetime

data = {"date": "", "debris": [], "spacecraft": []}
ref = "data/"


def randNum() -> float:
    return random.uniform(0.0, 100.0)


def create_spacecraft() -> dict:
    spacecraft = {}
    # spacecraft["id"] = uuid.uuid4().int & (1 >> 32)-1
    spacecraft["id"] = random.getrandbits(32)
    spacecraft["x_dis"] = randNum()
    spacecraft["y_dis"] = randNum()
    spacecraft["z_dis"] = randNum()
    spacecraft["x_vel"] = randNum()
    spacecraft["y_vel"] = randNum()
    spacecraft["z_vel"] = randNum()
    return spacecraft


def create_debris() -> dict:
    debris = {}
    # debris["id"] = uuid.uuid4().int & (1 >> 32)-1
    debris["id"] = random.getrandbits(32)
    debris["x_dis"] = randNum()
    debris["y_dis"] = randNum()
    debris["z_dis"] = randNum()
    debris["x_vel"] = randNum()
    debris["y_vel"] = randNum()
    debris["z_vel"] = randNum()
    return debris


def main(argv):

    output_file = ''
    debris_num = 0
    spacecraft_num = 0

    try:
        opts, args = getopt.getopt(argv, "hd:s:o:")
    except getopt.GetoptError:
        print('generate.py -h <help> -d <number of debris> -s <number of spacecraft> -o <output file name>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('generate.py -h <help> -d <number of debris> -s <number of spacecraft> -o <output file name>')
            sys.exit()
        elif opt in "-d":
            debris_num = int(arg)
        elif opt in "-s":
            spacecraft_num = int(arg)
        elif opt in "-o":
            output_file = arg

    data["date"] = datetime.datetime.utcnow().isoformat()

    for _ in range(debris_num):
        data["debris"].append(create_debris())

    for _ in range(spacecraft_num):
        data["spacecraft"].append(create_spacecraft())

    if output_file == '':
        with open(ref + "test_input.json", "w") as write_file:
            json.dump(data, write_file)
    else:
        with open(ref + output_file, "w") as write_file:
            json.dump(data, write_file)


if __name__ == '__main__':
    main(sys.argv[1:])
