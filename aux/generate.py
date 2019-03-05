
import json
import random
import uuid

data = {"Debris": [], "Spacecraft": []}


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


for i in range(100):
    data["Debris"].append(create_debris())
    data["Spacecraft"].append(create_spacecraft())


ref = "data/"

with open(ref + "test_input.json", "w") as write_file:
    json.dump(data, write_file)

