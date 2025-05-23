import toml

CONFIG_PATH = ""


def load_config(config_path:str) -> dict:
    global CONFIG_PATH
    CONFIG_PATH = config_path
    return toml.load(config_path)


def update_nested_toml(section_path:str, key:str, value: str | int | float | bool) -> None:
    try:
        with open(CONFIG_PATH, "r") as f:
            config = toml.load(f)
        sections = section_path.split(".")
        current = config

        for section in sections:
            current = current.setdefault(section, {})
        current[key] = value

        with open(CONFIG_PATH, "w") as f:
            toml.dump(config, f)

    except Exception as e:
        print(f"An error occurred: {e}")
