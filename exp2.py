
def adjacent_settings(inpset : dict[str,int|str]) -> list[dict[str,int|str]]:
    out : list[dict[str,int|str]] = list()
    # mirror iso
    s = inpset.copy()
    s["iso"] = 1 - s["iso"]
    out.append(s)
    # mirror sup
    s = inpset.copy()
    s["sup"] = 1 - s["sup"]
    out.append(s)
    if inpset["dec"]:
        pass
    else:
        
