from molli import Concurrent
from asyncio import create_subprocess_shell, sleep


async def proc(x, k=2, b=5):
    p = await create_subprocess_shell(f"sleep 0.5")
    await p.wait()
    # await sleep(0.25)
    return k * x + b


aa = Concurrent(range(100), update=3)

r = aa(proc)(k=2, b=5)
